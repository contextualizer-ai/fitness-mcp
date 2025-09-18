################################################################################
# fitness_mcp/main.py
# FastMCP wrapper for Agrobacterium mutant fitness data analysis
#
# This MCP analyzes fitness data from barcoded Agrobacterium mutants grown in
# mixed cultures across different conditions. Each row represents a gene knockout
# mutant, and fitness scores indicate:
# - NEGATIVE values: Gene knockout REDUCES fitness (gene is beneficial/essential for growth)
# - POSITIVE values: Gene knockout IMPROVES fitness (gene normally inhibits growth)
# - Values near 0: Gene knockout has minimal effect on fitness
################################################################################
import csv
import os
import random
from collections import Counter, defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple
from dataclasses import dataclass
from fastmcp import FastMCP
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
from fitness_mcp.data_processing import generate_significant_fitness_pairs


# METADATA REGISTRY SECTION
import threading
from functools import lru_cache


@dataclass
class GeneMetadata:
    """Metadata for a gene/locus (no fitness data)."""

    locus_id: str
    sys_name: str
    description: str


@dataclass
class ConditionMetadata:
    """Metadata for an experimental condition."""

    condition_id: str
    short_desc: str
    long_desc: str
    media: str
    temperature: str
    pH: str
    aerobic: str
    condition_1: str
    concentration_1: str
    units_1: str
    exp_group: str


@dataclass
class ModuleMetadata:
    """Metadata for a functional module."""

    module_id: int
    name: str
    category: str
    gene_count: int
    gene_list: List[str]


@dataclass
class FitnessEffect:
    """A pre-filtered significant fitness effect."""

    gene_id: str
    condition_id: str
    fitness_value: float


class MetadataRegistry:
    """Centralized registry for all gene, condition, and module metadata.

    Provides efficient, unified access to metadata while keeping fitness data separate.
    Eliminates redundancy between the old FitnessDataLoader, ModuleDataLoader, and PairsDataLoader.
    """

    def __init__(
        self,
        fitness_file: str = "data/fit_t.tab",
        exp_desc_file: str = "data/exp_organism_Agro.txt",
        modules_file: str = "data/RbTnSeq_modules_t1e-7.csv",
        module_meta_file: str = "data/module_meta.tsv",
        pairs_file: str = "data/fit_t_pairs_threshold_2_long.tab",
    ):
        """Initialize the metadata registry.

        Args:
            fitness_file: Path to fitness data file (for gene/condition extraction)
            exp_desc_file: Path to experimental conditions description file
            modules_file: Path to gene-module assignments file
            module_meta_file: Path to module metadata file
            pairs_file: Path to significant fitness pairs file
        """
        self.fitness_file = fitness_file
        self.exp_desc_file = exp_desc_file
        self.modules_file = modules_file
        self.module_meta_file = module_meta_file
        self.pairs_file = pairs_file

        # Core data structures
        self.genes: Dict[str, GeneMetadata] = {}
        self.conditions: Dict[str, ConditionMetadata] = {}
        self.modules: Dict[int, ModuleMetadata] = {}
        self.fitness_effects: List[FitnessEffect] = []

        # Index structures for efficient lookup
        self.gene_to_modules: Dict[str, List[int]] = {}
        self.module_to_genes: Dict[int, List[str]] = {}
        self.gene_to_conditions: Dict[str, List[str]] = {}
        self.condition_to_genes: Dict[str, List[str]] = {}

        # File tracking and thread safety
        self.loaded = False
        self._file_mtimes: Dict[str, float] = {}
        self._lock = threading.RLock()

    def _needs_reload(self) -> bool:
        """Check if any data file has been modified since last load."""
        try:
            files_to_check = [
                self.fitness_file,
                self.exp_desc_file,
                self.modules_file,
                self.module_meta_file,
                self.pairs_file,
            ]

            for filepath in files_to_check:
                if os.path.exists(filepath):
                    current_mtime = os.path.getmtime(filepath)
                    if (
                        filepath not in self._file_mtimes
                        or self._file_mtimes[filepath] != current_mtime
                    ):
                        return True

            return not self.loaded

        except OSError:
            return False

    def load_data(self) -> None:
        """Load all metadata from files with thread safety."""
        with self._lock:
            if not self._needs_reload():
                return

            # Clear existing data
            self._clear_data()

            # Load in dependency order
            self._load_gene_metadata()
            self._load_condition_metadata()
            self._load_module_metadata()
            self._load_fitness_effects()

            # Build index structures
            self._build_indexes()

            # Update file timestamps
            self._update_file_mtimes()

            self.loaded = True

    def _clear_data(self) -> None:
        """Clear all data structures."""
        self.genes.clear()
        self.conditions.clear()
        self.modules.clear()
        self.fitness_effects.clear()
        self.gene_to_modules.clear()
        self.module_to_genes.clear()
        self.gene_to_conditions.clear()
        self.condition_to_genes.clear()

    def _load_gene_metadata(self) -> None:
        """Load gene metadata from fitness file (first 3 columns only)."""
        if not os.path.exists(self.fitness_file):
            return

        with open(self.fitness_file, "r") as f:
            reader = csv.reader(f, delimiter="\t")

            # Skip header
            next(reader)

            # Read gene metadata (first 3 columns: locusId, sysName, desc)
            for row in reader:
                if len(row) < 3:
                    continue

                locus_id = row[0]
                sys_name = row[1]
                description = row[2]

                gene_meta = GeneMetadata(
                    locus_id=locus_id, sys_name=sys_name, description=description
                )

                # Index by both locus_id and sys_name
                self.genes[locus_id] = gene_meta
                if sys_name and sys_name != locus_id:
                    self.genes[sys_name] = gene_meta

    def _load_condition_metadata(self) -> None:
        """Load condition metadata from experimental description file."""
        if not os.path.exists(self.exp_desc_file):
            return

        with open(self.exp_desc_file, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader)

            for row in reader:
                if len(row) < len(header):
                    continue

                exp_data = dict(zip(header, row))
                exp_name = exp_data.get("expName", "")

                if exp_name:
                    condition_meta = ConditionMetadata(
                        condition_id=exp_name,
                        short_desc=exp_data.get("expDesc", ""),
                        long_desc=exp_data.get("expDescLong", ""),
                        media=exp_data.get("media", ""),
                        temperature=exp_data.get("temperature", ""),
                        pH=exp_data.get("pH", ""),
                        aerobic=exp_data.get("aerobic", ""),
                        condition_1=exp_data.get("condition_1", ""),
                        concentration_1=exp_data.get("concentration_1", ""),
                        units_1=exp_data.get("units_1", ""),
                        exp_group=exp_data.get("expGroup", ""),
                    )
                    self.conditions[exp_name] = condition_meta

    def _load_module_metadata(self) -> None:
        """Load module metadata and gene-module assignments."""
        # Load module metadata first
        if os.path.exists(self.module_meta_file):
            with open(self.module_meta_file, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    module_id = int(row["module"])
                    self.modules[module_id] = ModuleMetadata(
                        module_id=module_id,
                        name=row["name"],
                        category=row["category"],
                        gene_count=int(row["count"]),
                        gene_list=[],  # Will be populated from assignments file
                    )

        # Load gene-module assignments
        if os.path.exists(self.modules_file):
            with open(self.modules_file, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    module_id = int(row["module"])
                    locus_tag = row["locus_tag"]

                    # Add to module's gene list
                    if module_id in self.modules:
                        self.modules[module_id].gene_list.append(locus_tag)

    def _load_fitness_effects(self) -> None:
        """Load significant fitness effects from pairs file."""
        if not os.path.exists(self.pairs_file):
            # Generate if it doesn't exist
            self._generate_pairs_file()

        if not os.path.exists(self.pairs_file):
            return

        with open(self.pairs_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")

            for row in reader:
                effect = FitnessEffect(
                    gene_id=row["gene_id"],
                    condition_id=row["condition_id"],
                    fitness_value=float(row["value"]),
                )
                self.fitness_effects.append(effect)

    def _generate_pairs_file(self) -> None:
        """Generate pairs file if it doesn't exist."""
        # Extract threshold from filename
        filename = os.path.basename(self.pairs_file)
        if "threshold_" in filename:
            threshold_str = filename.split("threshold_")[1].split("_")[0]
            try:
                threshold = float(threshold_str)
            except ValueError:
                threshold = 2.0
        else:
            threshold = 2.0

        if os.path.exists(self.fitness_file):
            print(f"Generating pairs file with threshold {threshold}...")
            generate_significant_fitness_pairs(
                self.fitness_file, self.pairs_file, threshold
            )

    def _build_indexes(self) -> None:
        """Build efficient lookup indexes."""
        # Build gene-module indexes
        for module_id, module_meta in self.modules.items():
            for gene_id in module_meta.gene_list:
                if gene_id not in self.gene_to_modules:
                    self.gene_to_modules[gene_id] = []
                self.gene_to_modules[gene_id].append(module_id)

                if module_id not in self.module_to_genes:
                    self.module_to_genes[module_id] = []
                self.module_to_genes[module_id].append(gene_id)

        # Build gene-condition indexes from fitness effects
        for effect in self.fitness_effects:
            gene_id = effect.gene_id
            condition_id = effect.condition_id

            if gene_id not in self.gene_to_conditions:
                self.gene_to_conditions[gene_id] = []
            self.gene_to_conditions[gene_id].append(condition_id)

            if condition_id not in self.condition_to_genes:
                self.condition_to_genes[condition_id] = []
            self.condition_to_genes[condition_id].append(gene_id)

    def _update_file_mtimes(self) -> None:
        """Update file modification times."""
        files_to_check = [
            self.fitness_file,
            self.exp_desc_file,
            self.modules_file,
            self.module_meta_file,
            self.pairs_file,
        ]

        for filepath in files_to_check:
            if os.path.exists(filepath):
                self._file_mtimes[filepath] = os.path.getmtime(filepath)

    # Convenience methods for common access patterns

    def get_gene(self, gene_id: str) -> Optional[GeneMetadata]:
        """Get gene metadata by locus_id or sys_name."""
        self.load_data()
        return self.genes.get(gene_id)

    def get_condition(self, condition_id: str) -> Optional[ConditionMetadata]:
        """Get condition metadata by condition_id."""
        self.load_data()
        return self.conditions.get(condition_id)

    def get_module(self, module_id: int) -> Optional[ModuleMetadata]:
        """Get module metadata by module_id."""
        self.load_data()
        return self.modules.get(module_id)

    def search_genes(self, query: str, limit: int = 10) -> List[GeneMetadata]:
        """Search genes by name or description."""
        self.load_data()
        query = query.lower()
        matches = []
        seen_locus_ids = set()

        for gene_meta in self.genes.values():
            # Skip duplicates (since we index by both locus_id and sys_name)
            if gene_meta.locus_id in seen_locus_ids:
                continue

            if (
                query in gene_meta.locus_id.lower()
                or query in gene_meta.sys_name.lower()
                or query in gene_meta.description.lower()
            ):
                matches.append(gene_meta)
                seen_locus_ids.add(gene_meta.locus_id)

                if len(matches) >= limit:
                    break

        return matches

    def get_gene_modules(self, gene_id: str) -> List[ModuleMetadata]:
        """Get all modules containing a gene."""
        self.load_data()
        module_ids = self.gene_to_modules.get(gene_id, [])
        return [self.modules[mid] for mid in module_ids if mid in self.modules]

    def get_module_genes(self, module_id: int) -> List[GeneMetadata]:
        """Get all genes in a module."""
        self.load_data()
        gene_ids = self.module_to_genes.get(module_id, [])
        genes = []
        for gid in gene_ids:
            gene_meta = self.genes.get(gid)
            if gene_meta:
                genes.append(gene_meta)
        return genes

    def get_gene_fitness_effects(self, gene_id: str) -> List[FitnessEffect]:
        """Get all significant fitness effects for a gene."""
        self.load_data()
        return [effect for effect in self.fitness_effects if effect.gene_id == gene_id]

    def get_condition_fitness_effects(self, condition_id: str) -> List[FitnessEffect]:
        """Get all significant fitness effects for a condition."""
        self.load_data()
        return [
            effect
            for effect in self.fitness_effects
            if effect.condition_id == condition_id
        ]

    def search_modules(self, query: str, limit: int = 10) -> List[ModuleMetadata]:
        """Search modules by name or category."""
        self.load_data()
        query = query.lower()
        matches = []

        for module_meta in self.modules.values():
            if (
                query in module_meta.name.lower()
                or query in module_meta.category.lower()
            ):
                matches.append(module_meta)

                if len(matches) >= limit:
                    break

        return matches

    def get_all_conditions(self, condition_filter: Optional[str] = None) -> List[str]:
        """Get list of all condition IDs, optionally filtered."""
        self.load_data()
        condition_ids = list(self.conditions.keys())

        if condition_filter:
            filter_lower = condition_filter.lower()
            condition_ids = [
                cid for cid in condition_ids if filter_lower in cid.lower()
            ]

        return condition_ids

    def get_all_modules_list(self) -> List[Dict[str, Any]]:
        """Get list of all modules with basic info."""
        self.load_data()
        return [
            {
                "module_id": module_meta.module_id,
                "name": module_meta.name,
                "category": module_meta.category,
                "count": module_meta.gene_count,
            }
            for module_meta in self.modules.values()
        ]

    def interpret_fitness_score(self, fitness_score: float) -> Dict[str, Any]:
        """Interpret a fitness score in biological terms."""
        if fitness_score is None:
            return {
                "interpretation": "No data available",
                "effect": "unknown",
                "magnitude": "unknown",
            }

        abs_score = abs(fitness_score)

        if abs_score < 0.2:
            magnitude = "minimal"
        elif abs_score < 0.5:
            magnitude = "moderate"
        elif abs_score < 1.0:
            magnitude = "strong"
        else:
            magnitude = "very strong"

        if fitness_score < -0.1:
            effect = "gene_benefits_growth"
            interpretation = f"Gene knockout reduces fitness ({magnitude} effect). This gene is beneficial/essential for growth in this condition."
        elif fitness_score > 0.1:
            effect = "gene_inhibits_growth"
            interpretation = f"Gene knockout improves fitness ({magnitude} effect). This gene normally inhibits growth in this condition."
        else:
            effect = "neutral"
            interpretation = (
                "Gene knockout has minimal effect on fitness in this condition."
            )

        return {
            "interpretation": interpretation,
            "effect": effect,
            "magnitude": magnitude,
            "score": fitness_score,
        }


# Global metadata registry instance
metadata_registry = MetadataRegistry()


# Legacy data loaders (to be deprecated after migration)
# DATA LOADING SECTION


class ModuleDataLoader:
    def __init__(
        self,
        modules_file: str = "data/RbTnSeq_modules_t1e-7.csv",
        meta_file: str = "data/module_meta.tsv",
    ):
        """Initialize the module data loader.

        Args:
            modules_file: Path to the modules CSV file
            meta_file: Path to the module metadata TSV file
        """
        self.modules_file = modules_file
        self.meta_file = meta_file
        self.gene_to_modules: Dict[str, List[Dict[str, Any]]] = {}
        self.module_to_genes: Dict[int, List[Dict[str, Any]]] = {}
        self.module_meta: Dict[int, Dict[str, Any]] = {}
        self.loaded = False
        self._modules_mtime = -1.0
        self._meta_mtime = -1.0
        self._lock = threading.RLock()

    def _needs_reload(self) -> bool:
        """Check if either data file has been modified since last load."""
        try:
            modules_mtime = os.path.getmtime(self.modules_file)
            meta_mtime = os.path.getmtime(self.meta_file)
        except OSError:
            return False
        return (
            not self.loaded
            or modules_mtime != self._modules_mtime
            or meta_mtime != self._meta_mtime
        )

    def load_data(self) -> None:
        """Load the module data from files with thread safety."""
        with self._lock:
            if not self._needs_reload():
                return

            # Clear cache when reloading
            self._clear_caches()

            # Load module metadata first
            if not os.path.exists(self.meta_file):
                raise OSError(f"Module metadata file not found: {self.meta_file}")

        with open(self.meta_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                module_id = int(row["module"])
                self.module_meta[module_id] = {
                    "module_id": module_id,
                    "name": row["name"],
                    "category": row["category"],
                    "count": int(row["count"]),
                }

            # Load gene-module assignments
            if not os.path.exists(self.modules_file):
                raise OSError(f"Modules file not found: {self.modules_file}")

            with open(self.modules_file, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    module_id = int(row["module"])
                    locus_tag = row["locus_tag"]
                    gene_weight = (
                        float(row["gene_weight"]) if row["gene_weight"] else 0.0
                    )

                    gene_info = {
                        "locus_tag": locus_tag,
                        "module_id": module_id,
                        "gene_weight": gene_weight,
                        "product": row.get("product", ""),
                        "description": row.get("Description", ""),
                        "preferred_name": row.get("Preferred_name", ""),
                        "protein_names": row.get("Protein names", ""),
                        "go_terms": row.get("Gene Ontology (GO)", ""),
                    }

                    # Add module metadata to gene info
                    if module_id in self.module_meta:
                        gene_info.update(
                            {
                                "module_name": self.module_meta[module_id]["name"],
                                "module_category": self.module_meta[module_id][
                                    "category"
                                ],
                            }
                        )

                    # Index by gene
                    if locus_tag not in self.gene_to_modules:
                        self.gene_to_modules[locus_tag] = []
                    self.gene_to_modules[locus_tag].append(gene_info)

                    # Index by module
                    if module_id not in self.module_to_genes:
                        self.module_to_genes[module_id] = []
                    self.module_to_genes[module_id].append(gene_info)

            self.loaded = True
            self._modules_mtime = os.path.getmtime(self.modules_file)
            self._meta_mtime = os.path.getmtime(self.meta_file)

    def _clear_caches(self) -> None:
        """Clear LRU caches when data is reloaded."""
        self._cached_search_modules.cache_clear()

    @lru_cache(maxsize=128)
    def _cached_search_modules(
        self, query: str, limit: int, data_version: float
    ) -> List[Dict[str, Any]]:
        """Cached implementation of module search."""
        query = query.lower()
        matches = []

        for module_id, meta in self.module_meta.items():
            if query in meta["name"].lower() or query in meta["category"].lower():
                genes = self.module_to_genes.get(module_id, [])
                matches.append(
                    {"module": meta, "genes": genes, "gene_count": len(genes)}
                )

                if len(matches) >= limit:
                    break

        return matches

    def get_modules_for_gene(self, gene_id: str) -> List[Dict[str, Any]]:
        """Get module information for a specific gene.

        Args:
            gene_id: Gene locus tag

        Returns:
            List of modules containing this gene
        """
        self.load_data()
        return self.gene_to_modules.get(gene_id, [])

    def get_genes_in_module(self, module_id: int) -> Dict[str, Any]:
        """Get all genes in a specific module.

        Args:
            module_id: Module ID number

        Returns:
            Dictionary with module info and gene list
        """
        self.load_data()

        if module_id not in self.module_meta:
            return {"error": f"Module {module_id} not found"}

        genes = self.module_to_genes.get(module_id, [])

        return {
            "module": self.module_meta[module_id],
            "genes": genes,
            "gene_count": len(genes),
        }

    def search_modules_by_name(
        self, query: str, limit: int = 10
    ) -> List[Dict[str, Any]]:
        """Search modules by name or category.

        Args:
            query: Search term
            limit: Maximum results

        Returns:
            List of matching modules with gene lists
        """
        self.load_data()
        return self._cached_search_modules(query, limit, self._modules_mtime)

    def get_all_modules(self) -> List[Dict[str, Any]]:
        """Get list of all modules with basic info.

        Returns:
            List of all modules
        """
        self.load_data()
        return list(self.module_meta.values())


class FitnessDataLoader:
    def __init__(
        self,
        data_file: str = "data/fit_t.tab",
        exp_desc_file: str = "data/exp_organism_Agro.txt",
    ):
        """Initialize the fitness data loader.

        Args:
            data_file: Path to the fitness data file
            exp_desc_file: Path to the experimental conditions description file
        """
        self.data_file = data_file
        self.exp_desc_file = exp_desc_file
        self.genes: Dict[str, Dict[str, Any]] = {}
        self.conditions: List[str] = []
        self.condition_details: Dict[str, Dict[str, str]] = {}
        self.loaded = False
        self._mtime = -1.0
        self._exp_mtime = -1.0
        self._lock = threading.RLock()

    def _needs_reload(self) -> bool:
        """Check if either data file has been modified since last load."""
        try:
            data_mtime = os.path.getmtime(self.data_file)
            exp_mtime = os.path.getmtime(self.exp_desc_file)
        except OSError:
            return False
        return (
            not self.loaded or data_mtime != self._mtime or exp_mtime != self._exp_mtime
        )

    def load_data(self) -> None:
        """Load the fitness data from the tab-separated file with thread safety."""
        with self._lock:
            if not self._needs_reload():
                return

            if not os.path.exists(self.data_file):
                raise FileNotFoundError(
                    f"Fitness data file not found: {self.data_file}"
                )

            # Clear cache when reloading
            self._clear_caches()

            # Load experimental conditions descriptions first
            self._load_condition_details()

            with open(self.data_file, "r") as f:
                reader = csv.reader(f, delimiter="\t")

                # Read header (condition names)
                header = next(reader)
                self.conditions = header[3:]  # Skip locusId, sysName, desc columns

                # Read gene data
                for row in reader:
                    if len(row) < 4:  # Skip incomplete rows
                        continue

                    locus_id = row[0]
                    sys_name = row[1]
                    description = row[2]
                    fitness_values = [
                        float(val) if val and val != "NA" else None for val in row[3:]
                    ]

                    self.genes[locus_id] = {
                        "locusId": locus_id,
                        "sysName": sys_name,
                        "description": description,
                        "fitness_values": fitness_values,
                    }

                    # Also index by sysName if different from locusId
                    if sys_name and sys_name != locus_id:
                        self.genes[sys_name] = self.genes[locus_id]

            self.loaded = True
            self._mtime = os.path.getmtime(self.data_file)
            self._exp_mtime = os.path.getmtime(self.exp_desc_file)

    def _load_condition_details(self) -> None:
        """Load experimental condition details from the description file."""
        if not os.path.exists(self.exp_desc_file):
            return  # Skip if file doesn't exist

        with open(self.exp_desc_file, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader)  # Read header

            for row in reader:
                if len(row) < len(header):
                    continue

                # Create a dictionary mapping header to values
                exp_data = dict(zip(header, row))
                exp_name = exp_data.get("expName", "")

                if exp_name:
                    self.condition_details[exp_name] = {
                        "short_desc": exp_data.get("expDesc", ""),
                        "long_desc": exp_data.get("expDescLong", ""),
                        "media": exp_data.get("media", ""),
                        "temperature": exp_data.get("temperature", ""),
                        "pH": exp_data.get("pH", ""),
                        "aerobic": exp_data.get("aerobic", ""),
                        "condition_1": exp_data.get("condition_1", ""),
                        "concentration_1": exp_data.get("concentration_1", ""),
                        "units_1": exp_data.get("units_1", ""),
                        "exp_group": exp_data.get("expGroup", ""),
                    }

    def _clear_caches(self) -> None:
        """Clear LRU caches when data is reloaded."""
        self._cached_search_genes.cache_clear()
        self._cached_get_conditions.cache_clear()

    def get_gene_info(self, gene_id: str) -> Optional[Dict[str, Any]]:
        """Get basic information about a gene.

        Args:
            gene_id: Gene locus ID or system name

        Returns:
            Dictionary with gene information or None if not found
        """
        self.load_data()
        return self.genes.get(gene_id)

    def get_gene_fitness(
        self, gene_id: str, condition_filter: Optional[str] = None
    ) -> Dict[str, Any]:
        """Get fitness data for a gene knockout mutant across growth conditions.

        Args:
            gene_id: Gene locus ID or system name
            condition_filter: Optional filter to match condition names (case-insensitive)

        Returns:
            Dictionary with gene info and fitness scores. Fitness scores represent:
            - Negative values: Gene knockout reduces fitness (gene is beneficial/essential for growth)
            - Positive values: Gene knockout improves fitness (gene normally inhibits growth)
            - Values near 0: Gene knockout has minimal effect on fitness
        """
        self.load_data()

        gene_info = self.genes.get(gene_id)
        if not gene_info:
            return {"error": f"Gene {gene_id} not found"}

        # Filter conditions if requested
        if condition_filter:
            condition_filter = condition_filter.lower()
            matching_conditions = []
            matching_values = []

            for i, condition in enumerate(self.conditions):
                if condition_filter in condition.lower():
                    matching_conditions.append(condition)
                    if i < len(gene_info["fitness_values"]):
                        matching_values.append(gene_info["fitness_values"][i])
                    else:
                        matching_values.append(None)
        else:
            matching_conditions = self.conditions
            matching_values = gene_info["fitness_values"]

        # Create condition-value pairs with details
        fitness_data = []
        for condition, value in zip(matching_conditions, matching_values):
            condition_info = {"condition": condition, "fitness": value}

            # Add condition details if available
            if condition in self.condition_details:
                details = self.condition_details[condition]
                condition_info["description"] = details["short_desc"]
                condition_info["experimental_group"] = details["exp_group"]

            fitness_data.append(condition_info)

        return {
            "gene": {
                "locusId": gene_info["locusId"],
                "sysName": gene_info["sysName"],
                "description": gene_info["description"],
            },
            "fitness_data": fitness_data,
            "total_conditions": len(fitness_data),
        }

    @lru_cache(maxsize=256)
    def _cached_search_genes(
        self, query: str, limit: int, data_version: float
    ) -> List[Dict[str, Any]]:
        """Cached implementation of gene search."""
        query = query.lower()
        matches = []
        seen_locus_ids = set()

        for gene_id, gene_info in self.genes.items():
            # Skip duplicates (since we index by both locusId and sysName)
            if gene_info["locusId"] in seen_locus_ids:
                continue

            # Check if query matches locusId, sysName, or description
            if (
                query in gene_info["locusId"].lower()
                or query in gene_info["sysName"].lower()
                or query in gene_info["description"].lower()
            ):
                matches.append(
                    {
                        "locusId": gene_info["locusId"],
                        "sysName": gene_info["sysName"],
                        "description": gene_info["description"],
                    }
                )
                seen_locus_ids.add(gene_info["locusId"])

                if len(matches) >= limit:
                    break

        return matches

    def search_genes(self, query: str, limit: int = 10) -> List[Dict[str, Any]]:
        """Search for genes by name or description.

        Args:
            query: Search term to match against gene names or descriptions
            limit: Maximum number of results to return

        Returns:
            List of matching gene information
        """
        self.load_data()
        return self._cached_search_genes(query, limit, self._mtime)

    @lru_cache(maxsize=128)
    def _cached_get_conditions(
        self, condition_filter: Optional[str], data_version: float
    ) -> List[str]:
        """Cached implementation of get_conditions."""
        if condition_filter:
            condition_filter = condition_filter.lower()
            return [
                cond for cond in self.conditions if condition_filter in cond.lower()
            ]
        else:
            return self.conditions.copy()

    def get_conditions(self, condition_filter: Optional[str] = None) -> List[str]:
        """Get list of available growth conditions.

        Args:
            condition_filter: Optional filter to match condition names

        Returns:
            List of condition names
        """
        self.load_data()
        return self._cached_get_conditions(condition_filter, self._mtime)

    def interpret_fitness_score(self, fitness_score: float) -> Dict[str, Any]:
        """Interpret a fitness score in biological terms.

        Args:
            fitness_score: Numerical fitness score

        Returns:
            Dictionary with interpretation of the fitness effect
        """
        if fitness_score is None:
            return {
                "interpretation": "No data available",
                "effect": "unknown",
                "magnitude": "unknown",
            }

        abs_score = abs(fitness_score)

        if abs_score < 0.2:
            magnitude = "minimal"
        elif abs_score < 0.5:
            magnitude = "moderate"
        elif abs_score < 1.0:
            magnitude = "strong"
        else:
            magnitude = "very strong"

        if fitness_score < -0.1:
            effect = "gene_benefits_growth"
            interpretation = f"Gene knockout reduces fitness ({magnitude} effect). This gene is beneficial/essential for growth in this condition."
        elif fitness_score > 0.1:
            effect = "gene_inhibits_growth"
            interpretation = f"Gene knockout improves fitness ({magnitude} effect). This gene normally inhibits growth in this condition."
        else:
            effect = "neutral"
            interpretation = (
                "Gene knockout has minimal effect on fitness in this condition."
            )

        return {
            "interpretation": interpretation,
            "effect": effect,
            "magnitude": magnitude,
            "score": fitness_score,
        }

    def get_condition_details(self, condition_name: str) -> Dict[str, Any]:
        """Get detailed information about a specific growth condition.

        Args:
            condition_name: Name of the condition (e.g., 'set10IT004')

        Returns:
            Dictionary with detailed condition information
        """
        self.load_data()

        if condition_name not in self.condition_details:
            return {"error": f"Condition {condition_name} not found"}

        details = self.condition_details[condition_name]
        return {
            "condition_name": condition_name,
            "short_description": details["short_desc"],
            "long_description": details["long_desc"],
            "experimental_group": details["exp_group"],
            "growth_conditions": {
                "media": details["media"],
                "temperature": str(details["temperature"]) + "°C"
                if details["temperature"]
                else "",
                "pH": details["pH"],
                "aerobic": details["aerobic"],
            },
            "treatment": {
                "compound": details["condition_1"],
                "concentration": details["concentration_1"],
                "units": details["units_1"],
            }
            if details["condition_1"]
            else None,
        }


class PairsDataLoader:
    """
    Loads gene-condition pairs with significant fitness values from a tab-separated file.

    Only pairs where the absolute fitness value (|value|) exceeds a threshold of 2 are included.
    This threshold is used to identify gene knockouts or conditions that have a strong effect
    on fitness, filtering out minor or insignificant changes. The data file
    'fit_t_pairs_threshold_2_long.tab' is pre-filtered to include only these significant pairs.

    Attributes:
        data_file: Path to the pairs data file (default: "data/fit_t_pairs_threshold_2_long.tab")
        gene_to_conditions: Mapping from gene locus tag to list of conditions with significant fitness
        condition_to_genes: Mapping from condition to list of genes with significant fitness
        loaded: Whether the data has been loaded
        _mtime: Last modification time of the data file
        _lock: Threading lock for safe concurrent access
    """

    def __init__(self, data_file: str = "data/fit_t_pairs_threshold_2_long.tab"):
        """Initialize the pairs data loader.

        Args:
            data_file: Path to the pairs data file
        """
        self.data_file = data_file
        self.gene_to_conditions: Dict[str, List[Dict[str, Any]]] = {}
        self.condition_to_genes: Dict[str, List[Dict[str, Any]]] = {}
        self.loaded = False
        self._mtime = -1.0
        self._lock = threading.RLock()

    def _needs_reload(self) -> bool:
        """Check if the data file has been modified since last load."""
        try:
            mtime = os.path.getmtime(self.data_file)
        except OSError:
            return False
        return not self.loaded or mtime != self._mtime

    def load_data(self) -> None:
        """Load the pairs data from the tab-separated file."""
        with self._lock:
            if not self._needs_reload():
                return

            # Generate pairs file if it doesn't exist
            if not os.path.exists(self.data_file):
                self._generate_pairs_file()

            if not os.path.exists(self.data_file):
                raise FileNotFoundError(f"Pairs data file not found: {self.data_file}")

            # Clear existing data
            self.gene_to_conditions.clear()
            self.condition_to_genes.clear()

            with open(self.data_file, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")

                for row in reader:
                    gene_id = row["gene_id"]
                    condition_id = row["condition_id"]
                    value = float(row["value"])

                    # Index by gene
                    if gene_id not in self.gene_to_conditions:
                        self.gene_to_conditions[gene_id] = []
                    self.gene_to_conditions[gene_id].append(
                        {"condition": condition_id, "value": value}
                    )

                    # Index by condition
                    if condition_id not in self.condition_to_genes:
                        self.condition_to_genes[condition_id] = []
                    self.condition_to_genes[condition_id].append(
                        {"gene": gene_id, "value": value}
                    )

            self.loaded = True
            self._mtime = os.path.getmtime(self.data_file)

    def _generate_pairs_file(self) -> None:
        """Generate the pairs file from the main fitness data if it doesn't exist."""
        # Extract threshold from filename (e.g., "fit_t_pairs_threshold_2_long.tab")
        filename = os.path.basename(self.data_file)
        if "threshold_" in filename:
            threshold_str = filename.split("threshold_")[1].split("_")[0]
            try:
                threshold = float(threshold_str)
            except ValueError:
                threshold = 2.0  # Default threshold
        else:
            threshold = 2.0

        # Source fitness data file
        fitness_file = os.path.join(os.path.dirname(self.data_file), "fit_t.tab")

        if os.path.exists(fitness_file):
            print(f"Generating pairs file with threshold {threshold}...")
            generate_significant_fitness_pairs(fitness_file, self.data_file, threshold)
        else:
            raise FileNotFoundError(
                f"Source fitness data file not found: {fitness_file}"
            )

    def get_conditions_for_gene(self, gene_id: str) -> List[Dict[str, Any]]:
        """Get all conditions where a gene has significant fitness values (|value| > 2).

        Args:
            gene_id: Gene locus tag

        Returns:
            List of conditions with their fitness values
        """
        self.load_data()
        return self.gene_to_conditions.get(gene_id, [])

    def get_genes_for_condition(self, condition_id: str) -> List[Dict[str, Any]]:
        """Get all genes with significant fitness values (|value| > 2) for a condition.

        Args:
            condition_id: Condition identifier

        Returns:
            List of genes with their fitness values
        """
        self.load_data()
        return self.condition_to_genes.get(condition_id, [])
    
    def get_genes_with_all_conditions(self, conditions: Set[str]) -> Set[str]:
        """Find genes that have pairs with ALL the given conditions."""
        if not conditions:
            return set()
        
        # Start with genes from first condition, then intersect with others
        all_genes = None
        for cond in conditions:
            cond_genes = {g["gene"] for g in self.get_genes_for_condition(cond)}
            if all_genes is None:
                all_genes = cond_genes
            else:
                all_genes &= cond_genes
        
        return all_genes or set()
    
    def get_conditions_with_all_genes(self, genes: Set[str]) -> Set[str]:
        """Find conditions that have pairs with ALL the given genes."""
        if not genes:
            return set()
        
        # Start with conditions from first gene, then intersect with others
        all_conditions = None
        for gene in genes:
            gene_conds = {c["condition"] for c in self.get_conditions_for_gene(gene)}
            if all_conditions is None:
                all_conditions = gene_conds
            else:
                all_conditions &= gene_conds
        
        return all_conditions or set()
    
    def build_iterative_module_internal(self, start_gene: str, start_condition: str, max_size: int = 10) -> Dict[str, Any]:
        """
        Internal method to build a tight functional module starting from a gene-condition pair.
        
        Uses an iterative approach that alternates between:
        1. Adding genes that have significant fitness values in ALL current conditions
        2. Adding conditions where ALL current genes have significant fitness values
        
        Args:
            start_gene: Starting gene locus ID
            start_condition: Starting condition ID
            max_size: Maximum total elements (genes + conditions) in module
            
        Returns:
            Dict containing the built module or error information
        """
        self.load_data()
        
        # Verify starting pair exists
        if (start_gene, start_condition) not in {(g["gene"], start_condition) for g in self.get_genes_for_condition(start_condition)}:
            return {
                "error": f"No significant fitness value found for gene {start_gene} in condition {start_condition}",
                "reason": "Starting gene-condition pair not found in pairs data"
            }
        
        module_genes = {start_gene}
        module_conditions = {start_condition}
        
        iteration = 0
        max_iterations = max_size * 2
        history = []
        
        while iteration < max_iterations:
            iteration += 1
            
            if iteration % 2 == 1:
                # Odd: Add a gene that has pairs with ALL current conditions
                candidate_genes = self.get_genes_with_all_conditions(module_conditions)
                candidate_genes -= module_genes
                
                if not candidate_genes:
                    history.append(f"Iteration {iteration}: No new genes found with all conditions")
                    break
                
                # Choose gene with highest average absolute fitness across conditions
                best_gene = None
                best_score = -1
                
                for gene in candidate_genes:
                    gene_data = self.get_conditions_for_gene(gene)
                    relevant_values = [
                        abs(c["value"]) for c in gene_data 
                        if c["condition"] in module_conditions
                    ]
                    if relevant_values:
                        avg_score = sum(relevant_values) / len(relevant_values)
                        if avg_score > best_score:
                            best_score = avg_score
                            best_gene = gene
                
                if best_gene:
                    module_genes.add(best_gene)
                    history.append(f"Iteration {iteration}: Added gene {best_gene} (avg fitness: {best_score:.2f})")
            
            else:
                # Even: Add a condition that has pairs with ALL current genes
                candidate_conditions = self.get_conditions_with_all_genes(module_genes)
                candidate_conditions -= module_conditions
                
                if not candidate_conditions:
                    history.append(f"Iteration {iteration}: No new conditions found with all genes")
                    break
                
                # Choose condition with highest average absolute fitness across genes
                best_condition = None
                best_score = -1
                
                for condition in candidate_conditions:
                    cond_data = self.get_genes_for_condition(condition)
                    relevant_values = [
                        abs(g["value"]) for g in cond_data 
                        if g["gene"] in module_genes
                    ]
                    if relevant_values:
                        avg_score = sum(relevant_values) / len(relevant_values)
                        if avg_score > best_score:
                            best_score = avg_score
                            best_condition = condition
                
                if best_condition:
                    module_conditions.add(best_condition)
                    history.append(f"Iteration {iteration}: Added condition {best_condition} (avg fitness: {best_score:.2f})")
            
            # Stop if we've reached max size - this indicates runaway growth
            if len(module_genes) + len(module_conditions) >= max_size:
                history.append(f"Iteration {iteration}: Reached maximum size limit - runaway growth detected")
                return {
                    "error": f"Module building terminated due to runaway growth (>{max_size} elements)",
                    "reason": "Large modules indicate overly broad relationships rather than tight functional clusters",
                    "suggestion": "Try starting with a more specific gene-condition pair or reduce max_size parameter",
                    "start_gene": start_gene,
                    "start_condition": start_condition,
                    "partial_build": {
                        "genes_found": len(module_genes),
                        "conditions_found": len(module_conditions),
                        "iterations": iteration,
                        "history": history[-3:]  # Show last few steps
                    }
                }
        
        # Require at least 2 genes and 2 conditions for a meaningful module
        if len(module_genes) < 2 or len(module_conditions) < 2:
            return {
                "error": "Module too small to be meaningful",
                "reason": f"Found only {len(module_genes)} genes and {len(module_conditions)} conditions",
                "suggestion": "This gene-condition pair may be too isolated or specific",
                "start_gene": start_gene,
                "start_condition": start_condition,
                "build_info": {
                    "iterations": iteration,
                    "history": history
                }
            }
        
        # Calculate module density (how tight the relationships are)
        density = len(module_genes) * len(module_conditions) / (len(module_genes) + len(module_conditions))
        
        # Add module categorization
        categorization = self.categorize_module(sorted(module_genes), sorted(module_conditions))
        
        return {
            "start_gene": start_gene,
            "start_condition": start_condition,
            "module": {
                "genes": sorted(module_genes),
                "conditions": sorted(module_conditions),
                "num_genes": len(module_genes),
                "num_conditions": len(module_conditions),
                "total_possible_pairs": len(module_genes) * len(module_conditions),
                "density": round(density, 2)
            },
            "categorization": categorization,
            "build_info": {
                "iterations": iteration,
                "max_size_limit": max_size,
                "stopped_naturally": True,
                "history": history
            },
            "interpretation": f"Tight functional module with {len(module_genes)} genes and {len(module_conditions)} conditions (density: {density:.2f}). {categorization.get('summary_description', '')}"
        }
    
    def build_multiple_modules(self, num_modules: int = 5, max_size: int = 10, seed: Optional[int] = None) -> List[Dict[str, Any]]:
        """
        Build multiple modules using random starting points.
        
        Args:
            num_modules: Number of modules to build
            max_size: Maximum size per module
            seed: Random seed for reproducibility
            
        Returns:
            List of module dictionaries (successful and failed attempts)
        """
        self.load_data()
        
        if seed is not None:
            random.seed(seed)
        
        # Create all possible pairs for random selection
        all_pairs: Set[Tuple[str, str]] = set()
        for gene_id in self.gene_to_conditions:
            for condition_data in self.get_conditions_for_gene(gene_id):
                all_pairs.add((gene_id, condition_data["condition"]))
        
        modules = []
        used_starts: Set[Tuple[str, str]] = set()
        
        for i in range(num_modules):
            # Choose a random starting pair
            available_pairs = list(all_pairs - used_starts)
            if not available_pairs:
                break
            
            start_gene, start_condition = random.choice(available_pairs)
            used_starts.add((start_gene, start_condition))
            
            module = self.build_iterative_module_internal(start_gene, start_condition, max_size)
            modules.append(module)
        
        return modules
    
    def categorize_module(self, module_genes: List[str], module_conditions: List[str]) -> Dict[str, Any]:
        """
        Categorize a module based on its fitness characteristics.
        
        Args:
            module_genes: List of genes in the module
            module_conditions: List of conditions in the module
            
        Returns:
            Dict containing module categorization details
        """
        self.load_data()
        
        # Collect all fitness values for genes in the module
        fitness_values = []
        condition_fitness_map = defaultdict(list)
        gene_fitness_map = defaultdict(list)
        
        for gene in module_genes:
            gene_data = self.get_conditions_for_gene(gene)
            for condition_data in gene_data:
                if condition_data["condition"] in module_conditions:
                    value = condition_data["value"]
                    fitness_values.append(value)
                    condition_fitness_map[condition_data["condition"]].append(value)
                    gene_fitness_map[gene].append(value)
        
        if not fitness_values:
            return {
                "primary_category": "Unknown",
                "variability": "Unknown",
                "condition_diversity": "No Data",
                "mean_fitness": 0.0,
                "std_fitness": 0.0,
                "significant_conditions": 0,
                "top_conditions": [],
                "error": "No fitness data found for module"
            }
        
        # Calculate statistics
        if HAS_NUMPY:
            fitness_array = np.array(fitness_values)
            mean_fitness = float(np.mean(fitness_array))
            std_fitness = float(np.std(fitness_array))
        else:
            # Fallback calculations without numpy
            mean_fitness = sum(fitness_values) / len(fitness_values)
            std_fitness = (sum((x - mean_fitness) ** 2 for x in fitness_values) / len(fitness_values)) ** 0.5
        
        # Count conditions with significant fitness (|value| > 2.0)
        significant_conditions = len([c for c, values in condition_fitness_map.items() 
                                    if any(abs(v) > 2.0 for v in values)])
        
        # Determine primary category based on average fitness
        if mean_fitness > 2.0:
            primary_category = "Growth Enhancing"
        elif mean_fitness < -2.0:
            primary_category = "Growth Inhibiting"
        else:
            primary_category = "Neutral/Adaptive"
        
        # Subcategories based on variability
        if std_fitness > 3.0:
            variability = "High Variability"
        elif std_fitness > 1.5:
            variability = "Moderate Variability"
        else:
            variability = "Low Variability"
        
        # Condition diversity assessment
        if significant_conditions > 50:
            condition_diversity = "Broad Condition Response"
        elif significant_conditions > 20:
            condition_diversity = "Moderate Condition Response"
        else:
            condition_diversity = "Narrow Condition Response"
        
        # Top conditions by average absolute fitness
        condition_avg_fitness = {}
        for condition, values in condition_fitness_map.items():
            if values:
                condition_avg_fitness[condition] = sum(abs(v) for v in values) / len(values)
        
        top_conditions = sorted(condition_avg_fitness.items(), 
                              key=lambda x: x[1], reverse=True)[:5]
        top_conditions = [cond for cond, _ in top_conditions]
        
        # Generate summary description
        summary = f"{primary_category} module with {variability.lower()} and {condition_diversity.lower()}"
        
        return {
            "primary_category": primary_category,
            "variability": variability,
            "condition_diversity": condition_diversity,
            "mean_fitness": round(mean_fitness, 3),
            "std_fitness": round(std_fitness, 3),
            "significant_conditions": significant_conditions,
            "top_conditions": top_conditions,
            "summary_description": summary,
            "module_size": {
                "num_genes": len(module_genes),
                "num_conditions": len(module_conditions),
                "total_pairs": len(fitness_values)
            }
        }
    
    def generate_module_summary(self, modules: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Generate a summary of categorized modules.
        
        Args:
            modules: List of modules with categorization data
            
        Returns:
            Dict containing summary statistics across all modules
        """
        if not modules:
            return {"error": "No modules provided for summary"}
        
        successful_modules = [m for m in modules if "error" not in m and "module" in m]
        
        if not successful_modules:
            return {"error": "No successful modules found for summary"}
        
        # Categorize modules by primary category
        category_counts: Counter[str] = Counter()
        variability_counts: Counter[str] = Counter()
        diversity_counts: Counter[str] = Counter()
        
        total_genes = set()
        total_conditions = set()
        avg_densities = []
        
        for module in successful_modules:
            if "categorization" in module:
                cat = module["categorization"]
                category_counts[cat["primary_category"]] += 1
                variability_counts[cat["variability"]] += 1
                diversity_counts[cat["condition_diversity"]] += 1
            
            if "module" in module:
                mod = module["module"]
                total_genes.update(mod["genes"])
                total_conditions.update(mod["conditions"])
                avg_densities.append(mod["density"])
        
        avg_density = sum(avg_densities) / len(avg_densities) if avg_densities else 0
        
        return {
            "summary": {
                "total_modules": len(successful_modules),
                "unique_genes": len(total_genes),
                "unique_conditions": len(total_conditions),
                "average_module_density": round(avg_density, 3)
            },
            "categorization_summary": {
                "primary_categories": dict(category_counts),
                "variability_patterns": dict(variability_counts),
                "condition_diversity": dict(diversity_counts)
            },
            "interpretation": f"Analyzed {len(successful_modules)} modules covering {len(total_genes)} genes and {len(total_conditions)} conditions. "
                           f"Most common category: {category_counts.most_common(1)[0][0] if category_counts else 'None'}"
        }


# Global data loader instances
fitness_loader = FitnessDataLoader()
module_loader = ModuleDataLoader()
pairs_loader = PairsDataLoader()


# MCP TOOL SECTION
def get_gene_info(gene_id: str) -> Dict[str, Any]:
    """
    Get basic information about a gene.

    Args:
        gene_id: Gene locus ID (e.g., 'Atu0001') or system name

    Returns:
        Dict containing gene information including locusId, sysName, and description
    """
    gene_meta = metadata_registry.get_gene(gene_id)
    if gene_meta:
        return {
            "locusId": gene_meta.locus_id,
            "sysName": gene_meta.sys_name,
            "description": gene_meta.description,
        }
    else:
        return {"error": f"Gene {gene_id} not found"}


def get_gene_fitness(
    gene_id: str, condition_filter: Optional[str] = None
) -> Dict[str, Any]:
    """
    Get fitness data for a gene across different growth conditions.

    Args:
        gene_id: Gene locus ID (e.g., 'Atu0001') or system name
        condition_filter: Optional filter to match specific conditions (e.g., 'glucose', 'stress', 'antibiotic')

    Returns:
        Dict containing gene info and fitness data for matching conditions
    """
    # Note: Still using legacy loader for now as it contains fitness values
    # This will be refactored once fitness data is properly separated
    return fitness_loader.get_gene_fitness(gene_id, condition_filter)


def search_genes(query: str, limit: int = 3) -> List[Dict[str, Any]]:
    """
    Search for genes by name or description.

    Args:
        query: Search term to match against gene names or descriptions
        limit: Maximum number of results to return (default: 3)

    Returns:
        List of dictionaries containing matching gene information
    """
    gene_metas = metadata_registry.search_genes(query, limit)
    return [
        {
            "locusId": gene_meta.locus_id,
            "sysName": gene_meta.sys_name,
            "description": gene_meta.description,
        }
        for gene_meta in gene_metas
    ]


def get_growth_conditions(condition_filter: Optional[str] = None) -> List[str]:
    """
    Get list of available growth conditions.

    Args:
        condition_filter: Optional filter to match condition names (e.g., 'LB', 'glucose', 'metal')

    Returns:
        List of condition names
    """
    return metadata_registry.get_all_conditions(condition_filter)


def get_condition_details(condition_name: str) -> Dict[str, Any]:
    """
    Get detailed information about a specific growth condition.

    Args:
        condition_name: Name of the condition (e.g., 'set10IT004', 'set10IT020')

    Returns:
        Dict with detailed condition information including description, media, temperature, pH, and treatment details
    """
    condition_meta = metadata_registry.get_condition(condition_name)
    if not condition_meta:
        return {"error": f"Condition {condition_name} not found"}

    return {
        "condition_name": condition_meta.condition_id,
        "short_description": condition_meta.short_desc,
        "long_description": condition_meta.long_desc,
        "experimental_group": condition_meta.exp_group,
        "growth_conditions": {
            "media": condition_meta.media,
            "temperature": str(condition_meta.temperature) + "°C"
            if condition_meta.temperature
            else "",
            "pH": condition_meta.pH,
            "aerobic": condition_meta.aerobic,
        },
        "treatment": {
            "compound": condition_meta.condition_1,
            "concentration": condition_meta.concentration_1,
            "units": condition_meta.units_1,
        }
        if condition_meta.condition_1
        else None,
    }


def interpret_fitness_score(fitness_score: float) -> Dict[str, Any]:
    """
    Interpret a fitness score in biological terms.

    Args:
        fitness_score: Numerical fitness score from gene knockout experiment

    Returns:
        Dict with biological interpretation of the fitness effect:
        - Negative scores: Gene knockout reduces fitness (gene is beneficial/essential for growth)
        - Positive scores: Gene knockout improves fitness (gene normally inhibits growth)
        - Near zero: Gene knockout has minimal effect
    """
    return metadata_registry.interpret_fitness_score(fitness_score)


def find_essential_genes(
    condition_filter: Optional[str] = None,
    min_fitness_threshold: float = 0.5,
    limit: int = 5,
) -> List[Dict[str, Any]]:
    """
    Find genes that appear essential (high positive fitness scores when knocked out).

    Args:
        condition_filter: Optional filter to match specific conditions
        min_fitness_threshold: Minimum fitness score to consider essential (default: 0.5)
        limit: Maximum number of genes to return

    Returns:
        List of genes with high positive fitness scores indicating essentiality
    """
    # Note: Still using legacy fitness loader as it provides access to full fitness data
    # This function requires complete fitness matrices, not just significant effects
    fitness_loader.load_data()

    essential_genes = []
    seen_genes = set()

    for gene_id, gene_info in fitness_loader.genes.items():
        if gene_info["locusId"] in seen_genes:
            continue
        seen_genes.add(gene_info["locusId"])

        # Get fitness data for this gene
        fitness_data = fitness_loader.get_gene_fitness(gene_id, condition_filter)
        if "error" in fitness_data:
            continue

        # Find conditions where gene appears essential
        essential_conditions = []
        for condition_data in fitness_data["fitness_data"]:
            if (
                condition_data["fitness"] is not None
                and condition_data["fitness"] >= min_fitness_threshold
            ):
                essential_conditions.append(
                    {
                        "condition": condition_data["condition"],
                        "fitness_score": condition_data["fitness"],
                        "description": condition_data.get("description", ""),
                        "interpretation": fitness_loader.interpret_fitness_score(
                            condition_data["fitness"]
                        ),
                    }
                )

        # Sort by fitness score and limit conditions per gene
        essential_conditions.sort(key=lambda x: x["fitness_score"], reverse=True)
        essential_conditions = essential_conditions[:3]  # Max 3 conditions per gene

        if essential_conditions:
            essential_genes.append(
                {
                    "gene": fitness_data["gene"],
                    "essential_in_conditions": essential_conditions,
                    "num_essential_conditions": len(essential_conditions),
                }
            )

        if len(essential_genes) >= limit:
            break

    # Sort by number of conditions where gene appears essential
    essential_genes.sort(key=lambda x: x["num_essential_conditions"], reverse=True)
    return essential_genes


def find_growth_inhibitor_genes(
    condition_filter: Optional[str] = None,
    max_fitness_threshold: float = -0.5,
    limit: int = 5,
) -> List[Dict[str, Any]]:
    """
    Find genes that inhibit growth (negative fitness scores when knocked out).

    Args:
        condition_filter: Optional filter to match specific conditions
        max_fitness_threshold: Maximum fitness score to consider inhibitory (default: -0.5)
        limit: Maximum number of genes to return

    Returns:
        List of genes with negative fitness scores indicating they normally inhibit growth
    """
    # Note: Still using legacy fitness loader as it provides access to full fitness data
    # This function requires complete fitness matrices, not just significant effects
    fitness_loader.load_data()

    inhibitor_genes = []
    seen_genes = set()

    for gene_id, gene_info in fitness_loader.genes.items():
        if gene_info["locusId"] in seen_genes:
            continue
        seen_genes.add(gene_info["locusId"])

        # Get fitness data for this gene
        fitness_data = fitness_loader.get_gene_fitness(gene_id, condition_filter)
        if "error" in fitness_data:
            continue

        # Find conditions where gene inhibits growth
        inhibitor_conditions = []
        for condition_data in fitness_data["fitness_data"]:
            if (
                condition_data["fitness"] is not None
                and condition_data["fitness"] <= max_fitness_threshold
            ):
                inhibitor_conditions.append(
                    {
                        "condition": condition_data["condition"],
                        "fitness_score": condition_data["fitness"],
                        "description": condition_data.get("description", ""),
                        "interpretation": fitness_loader.interpret_fitness_score(
                            condition_data["fitness"]
                        ),
                    }
                )

        # Sort by fitness score and limit conditions per gene
        inhibitor_conditions.sort(key=lambda x: x["fitness_score"])
        inhibitor_conditions = inhibitor_conditions[:3]  # Max 3 conditions per gene

        if inhibitor_conditions:
            inhibitor_genes.append(
                {
                    "gene": fitness_data["gene"],
                    "inhibits_growth_in_conditions": inhibitor_conditions,
                    "num_inhibitory_conditions": len(inhibitor_conditions),
                }
            )

        if len(inhibitor_genes) >= limit:
            break

    # Sort by number of conditions where gene inhibits growth
    inhibitor_genes.sort(key=lambda x: x["num_inhibitory_conditions"], reverse=True)
    return inhibitor_genes


def analyze_gene_fitness(
    gene_id: str,
    min_fitness: Optional[float] = None,
    max_fitness: Optional[float] = None,
    limit: Optional[int] = 10,
) -> Dict[str, Any]:
    """
    Analyze fitness effects for a gene knockout mutant across conditions.

    Args:
        gene_id: Gene locus ID or system name
        min_fitness: Minimum fitness value to include (use negative values to find essential genes)
        max_fitness: Maximum fitness value to include (use positive values to find growth-inhibiting genes)
        limit: Maximum number of conditions to return per category (default: 10, sorted by absolute fitness)

    Returns:
        Dict with gene info and categorized fitness data. Categories indicate:
        - conditions_where_gene_is_essential: Negative fitness scores (gene knockout reduces fitness, gene is ESSENTIAL)
        - conditions_where_gene_inhibits_growth: Positive fitness scores (gene knockout improves fitness, gene INHIBITS growth)
        - neutral_conditions: Fitness scores near zero (gene knockout has minimal effect)
    """
    # Note: Still using legacy fitness loader as it provides access to full fitness data
    # This function requires complete fitness matrices, not just significant effects
    fitness_data = fitness_loader.get_gene_fitness(gene_id)

    if "error" in fitness_data:
        return fitness_data

    # Categorize conditions by fitness effect
    essential_conditions = []  # Negative fitness: gene knockout reduces fitness (gene is ESSENTIAL)
    inhibitory_conditions = []  # Positive fitness: gene knockout improves fitness (gene INHIBITS growth)
    neutral = []  # Near-zero fitness: gene knockout has minimal effect

    for item in fitness_data["fitness_data"]:
        fitness_val = item["fitness"]
        if fitness_val is None:
            continue

        # Apply filters if specified
        if min_fitness is not None and fitness_val < min_fitness:
            continue
        if max_fitness is not None and fitness_val > max_fitness:
            continue

        if (
            fitness_val < -0.5
        ):  # Gene knockout reduces fitness (gene is ESSENTIAL for this condition)
            essential_conditions.append(item)
        elif (
            fitness_val > 0.5
        ):  # Gene knockout improves fitness (gene INHIBITS growth in this condition)
            inhibitory_conditions.append(item)
        else:  # Neutral effect
            neutral.append(item)

    # Sort each category by fitness value and apply limit
    # Sort essential_conditions by most negative first (strongest essential effect)
    essential_conditions.sort(key=lambda x: x["fitness"])
    if limit:
        essential_conditions = essential_conditions[:limit]

    # Sort inhibitory_conditions by most positive first (strongest inhibitory effect)
    inhibitory_conditions.sort(key=lambda x: x["fitness"], reverse=True)
    if limit:
        inhibitory_conditions = inhibitory_conditions[:limit]

    # Sort neutral by absolute value (closest to zero first)
    neutral.sort(key=lambda x: abs(x["fitness"]))
    if limit:
        neutral = neutral[:limit]

    return {
        "gene": fitness_data["gene"],
        "analysis": {
            "conditions_where_gene_is_essential": essential_conditions,
            "conditions_where_gene_inhibits_growth": inhibitory_conditions,
            "neutral_conditions": neutral,
            "summary": {
                "total_conditions_tested": len(
                    [
                        x
                        for x in fitness_data["fitness_data"]
                        if x["fitness"] is not None
                    ]
                ),
                "essential_count": len(essential_conditions),
                "inhibitory_count": len(inhibitory_conditions),
                "neutral_count": len(neutral),
                "limit_applied": limit,
            },
        },
    }


def get_gene_modules(gene_id: str) -> Dict[str, Any]:
    """
    Get module information for a specific gene/locus.

    Args:
        gene_id: Gene locus tag (e.g., 'Atu0001')

    Returns:
        Dict containing modules that include this gene
    """
    module_metas = metadata_registry.get_gene_modules(gene_id)

    if not module_metas:
        return {"error": f"No modules found for gene {gene_id}"}

    modules = [
        {
            "locus_tag": gene_id,
            "module_id": module_meta.module_id,
            "gene_weight": 1.0,  # Default weight - could be enhanced later
            "product": "",  # Could be populated from gene description
            "module_name": module_meta.name,
            "module_category": module_meta.category,
        }
        for module_meta in module_metas
    ]

    return {"gene_id": gene_id, "modules": modules, "module_count": len(modules)}


def get_module_genes(module_id: int) -> Dict[str, Any]:
    """
    Get all genes in a specific module.

    Args:
        module_id: Module ID number

    Returns:
        Dict containing module info and all genes in the module
    """
    module_meta = metadata_registry.get_module(module_id)
    if not module_meta:
        return {"error": f"Module {module_id} not found"}

    gene_metas = metadata_registry.get_module_genes(module_id)

    genes = [
        {
            "locus_tag": gene_meta.locus_id,
            "module_id": module_id,
            "gene_weight": 1.0,  # Default weight
            "product": gene_meta.description,
        }
        for gene_meta in gene_metas
    ]

    return {
        "module": {
            "module_id": module_meta.module_id,
            "name": module_meta.name,
            "category": module_meta.category,
            "count": module_meta.gene_count,
        },
        "genes": genes,
        "gene_count": len(genes),
    }


def search_modules(query: str, limit: int = 3) -> List[Dict[str, Any]]:
    """
    Search for modules by name or category.

    Args:
        query: Search term to match against module names or categories
        limit: Maximum number of results to return (default: 3)

    Returns:
        List of matching modules with their genes
    """
    module_metas = metadata_registry.search_modules(query, limit)

    results = []
    for module_meta in module_metas:
        gene_metas = metadata_registry.get_module_genes(module_meta.module_id)

        genes = [
            {
                "locus_tag": gene_meta.locus_id,
                "module_id": module_meta.module_id,
                "gene_weight": 1.0,
                "product": gene_meta.description,
            }
            for gene_meta in gene_metas
        ]

        results.append(
            {
                "module": {
                    "module_id": module_meta.module_id,
                    "name": module_meta.name,
                    "category": module_meta.category,
                    "count": module_meta.gene_count,
                },
                "genes": genes,
                "gene_count": len(genes),
            }
        )

    return results


def get_all_modules() -> List[Dict[str, Any]]:
    """
    Get list of all available modules.

    Returns:
        List of all modules with basic information
    """
    return metadata_registry.get_all_modules_list()


def get_conditions_for_gene(gene_id: str) -> Dict[str, Any]:
    """
    Get all conditions where a gene has significant fitness values (|value| > 2).

    Args:
        gene_id: Gene locus ID (e.g., 'Atu0001')

    Returns:
        Dict containing gene_id and list of conditions with their fitness values
    """
    fitness_effects = metadata_registry.get_gene_fitness_effects(gene_id)

    if not fitness_effects:
        return {"error": f"No significant conditions found for gene {gene_id}"}

    # Convert to expected format and sort by absolute value of fitness score
    conditions = [
        {"condition": effect.condition_id, "value": effect.fitness_value}
        for effect in fitness_effects
    ]
    conditions.sort(
        key=lambda x: abs(x["value"]) if isinstance(x["value"], (int, float)) else 0,
        reverse=True,
    )

    return {
        "gene_id": gene_id,
        "conditions": conditions,
        "total_conditions": len(conditions),
        "interpretation": "These are conditions where the gene knockout has significant fitness effects (|value| > 2)",
    }


def get_genes_for_condition(condition_id: str) -> Dict[str, Any]:
    """
    Get all genes with significant fitness values (|value| > 2) for a condition.

    Args:
        condition_id: Condition identifier (e.g., 'set10IT004 D,L-Malic Acid (C)')

    Returns:
        Dict containing condition_id and list of genes with their fitness values
    """
    fitness_effects = metadata_registry.get_condition_fitness_effects(condition_id)

    if not fitness_effects:
        return {
            "error": f"No genes with significant fitness values found for condition {condition_id}"
        }

    # Convert to expected format and sort by absolute value of fitness score
    genes = [
        {"gene": effect.gene_id, "value": effect.fitness_value}
        for effect in fitness_effects
    ]
    genes.sort(
        key=lambda x: abs(x["value"]) if isinstance(x["value"], (int, float)) else 0,
        reverse=True,
    )

    return {
        "condition_id": condition_id,
        "genes": genes,
        "total_genes": len(genes),
        "interpretation": "These are genes where knockout has significant fitness effects (|value| > 2) in this condition",
    }


def expand_gene_condition_network(gene_id: str, condition_id: str) -> Dict[str, Any]:
    """
    Perform a two-hop expansion from a gene-condition pair to find related genes and conditions.

    Starting from a specific gene-condition pair, this function:
    1. Finds all conditions where the gene has significant fitness values
    2. Finds all genes that have significant fitness values in the condition
    3. Expands to find all conditions for the gene set and all genes for the condition set

    Args:
        gene_id: Gene locus ID (e.g., 'Atu0001')
        condition_id: Condition identifier (e.g., 'set10IT004 D,L-Malic Acid (C)')

    Returns:
        Dict containing the expanded network of related genes and conditions
    """
    metadata_registry.load_data()

    # Check if the gene-condition pair exists (use index for quick lookup)
    conditions_for_gene = set(metadata_registry.gene_to_conditions.get(gene_id, []))

    if condition_id not in conditions_for_gene:
        return {
            "error": f"No significant fitness value found for gene {gene_id} in condition {condition_id}"
        }

    # Get the specific fitness value for this pair (need to scan once)
    gene_effects = metadata_registry.get_gene_fitness_effects(gene_id)
    gene_conditions = {
        effect.condition_id: effect.fitness_value for effect in gene_effects
    }

    # Step 1: conditions_for_gene already set above from index

    # Step 2: Get all genes for this condition (first hop from condition - use index)
    genes_for_condition = set(
        metadata_registry.condition_to_genes.get(condition_id, [])
    )

    # Step 3: Expand - get all genes for the condition set (using index for efficiency)
    all_genes: set[str] = set()
    for cond in conditions_for_gene:
        # Use index structure instead of linear scan
        genes_in_cond = metadata_registry.condition_to_genes.get(cond, [])
        all_genes.update(genes_in_cond)

    # Step 4: Expand - get all conditions for the gene set (using index for efficiency)
    all_conditions: set[str] = set()
    for gene in genes_for_condition:
        # Use index structure instead of linear scan
        conditions_for_gene_item = metadata_registry.gene_to_conditions.get(gene, [])
        all_conditions.update(conditions_for_gene_item)

    # Get the specific fitness value for the query pair
    query_value = gene_conditions[condition_id]

    return {
        "query": {
            "gene_id": gene_id,
            "condition_id": condition_id,
            "fitness_value": query_value,
        },
        "first_hop": {
            "conditions_for_query_gene": sorted(conditions_for_gene),
            "genes_for_query_condition": sorted(genes_for_condition),
            "num_conditions": len(conditions_for_gene),
            "num_genes": len(genes_for_condition),
        },
        "second_hop": {
            "all_genes_in_network": sorted(all_genes),
            "all_conditions_in_network": sorted(all_conditions),
            "num_total_genes": len(all_genes),
            "num_total_conditions": len(all_conditions),
        },
        "network_size": {
            "gene_expansion_factor": len(all_genes) / len(genes_for_condition)
            if genes_for_condition
            else 0,
            "condition_expansion_factor": len(all_conditions) / len(conditions_for_gene)
            if conditions_for_gene
            else 0,
        },
        "interpretation": "This network shows genes and conditions related through significant fitness effects",
    }


def build_iterative_module(gene_id: str, condition_id: str, max_size: int = 10) -> Dict[str, Any]:
    """
    Build a tight functional module starting from a gene-condition pair.
    
    Uses an iterative approach that alternates between:
    1. Adding genes that have significant fitness values in ALL current conditions
    2. Adding conditions where ALL current genes have significant fitness values
    
    This creates coherent modules rather than large submatrices.
    
    Args:
        gene_id: Starting gene locus ID (e.g., 'Atu0001')
        condition_id: Starting condition ID
        max_size: Maximum total elements (genes + conditions) in module
        
    Returns:
        Dict containing the built module with genes, conditions, and build history
    """
    return pairs_loader.build_iterative_module_internal(gene_id, condition_id, max_size)


def discover_functional_modules(num_modules: int = 10, max_size: int = 10, seed: Optional[int] = None) -> Dict[str, Any]:
    """
    Discover multiple functional modules by trying random gene-condition starting points.
    
    This function attempts to build multiple tight functional modules using the iterative
    approach. It provides a summary of successful vs failed attempts and returns only
    the modules that represent tight, biologically meaningful relationships.
    
    Args:
        num_modules: Number of module building attempts (default: 10)
        max_size: Maximum size per module to prevent runaway growth (default: 10)
        seed: Random seed for reproducibility (default: None)
        
    Returns:
        Dict containing successful modules, failed attempts summary, and statistics
    """
    modules = pairs_loader.build_multiple_modules(num_modules, max_size, seed)
    
    successful_modules = [m for m in modules if "error" not in m]
    failed_modules = [m for m in modules if "error" in m]
    
    # Categorize failure reasons
    runaway_growth = len([m for m in failed_modules if "runaway growth" in m.get("error", "")])
    too_small = len([m for m in failed_modules if "too small" in m.get("error", "")])
    not_found = len([m for m in failed_modules if "not found" in m.get("error", "")])
    
    # Calculate statistics for successful modules
    if successful_modules:
        avg_genes = sum(m["module"]["num_genes"] for m in successful_modules) / len(successful_modules)
        avg_conditions = sum(m["module"]["num_conditions"] for m in successful_modules) / len(successful_modules)
        avg_density = sum(m["module"]["density"] for m in successful_modules) / len(successful_modules)
        avg_iterations = sum(m["build_info"]["iterations"] for m in successful_modules) / len(successful_modules)
    else:
        avg_genes = avg_conditions = avg_density = avg_iterations = 0
    
    # Generate categorization summary
    categorization_summary = pairs_loader.generate_module_summary(successful_modules)
    
    return {
        "summary": {
            "total_attempts": num_modules,
            "successful_modules": len(successful_modules),
            "failed_attempts": len(failed_modules),
            "success_rate": round(len(successful_modules) / num_modules, 2) if num_modules > 0 else 0
        },
        "failure_analysis": {
            "runaway_growth": runaway_growth,
            "too_small": too_small,
            "pair_not_found": not_found
        },
        "module_statistics": {
            "average_genes_per_module": round(avg_genes, 1),
            "average_conditions_per_module": round(avg_conditions, 1),
            "average_density": round(avg_density, 2),
            "average_iterations": round(avg_iterations, 1)
        },
        "categorization_analysis": categorization_summary,
        "successful_modules": successful_modules,
        "parameters": {
            "max_size_limit": max_size,
            "random_seed": seed
        },
        "interpretation": f"Found {len(successful_modules)} tight functional modules out of {num_modules} attempts. "
                         f"Modules represent coherent gene-condition relationships with average density {avg_density:.2f}. "
                         f"{categorization_summary.get('interpretation', '')}"
    }


# MAIN SECTION
# Create the FastMCP instance
mcp = FastMCP("fitness_mcp")

# Register all tools
mcp.tool(get_gene_info)
mcp.tool(get_gene_fitness)
mcp.tool(search_genes)
mcp.tool(get_growth_conditions)
mcp.tool(get_condition_details)
mcp.tool(interpret_fitness_score)
mcp.tool(find_essential_genes)
mcp.tool(find_growth_inhibitor_genes)
mcp.tool(analyze_gene_fitness)
mcp.tool(get_gene_modules)
mcp.tool(get_module_genes)
mcp.tool(search_modules)
mcp.tool(get_all_modules)
mcp.tool(get_conditions_for_gene)
mcp.tool(get_genes_for_condition)
mcp.tool(expand_gene_condition_network)
mcp.tool(build_iterative_module)
mcp.tool(discover_functional_modules)


def main() -> None:
    """Main entry point for the application."""
    mcp.run()


if __name__ == "__main__":
    main()
