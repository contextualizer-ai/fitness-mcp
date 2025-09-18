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
from operator import itemgetter
from typing import Any, Dict, List, Optional
from dataclasses import dataclass
from fastmcp import FastMCP
from fitness_mcp.data_processing import generate_significant_fitness_pairs


# CONFIGURATION SECTION
class FitnessConfig:
    """Configuration for fitness analysis thresholds and limits.

    Centralizes all configurable parameters for POSIX-style tool consistency.

    Note: __slots__ is used here for Python memory optimization (not LinkML schema).
    It restricts instance attributes to only those listed, reducing memory overhead
    and improving performance for frequently instantiated objects.
    """

    __slots__ = (
        "significance_threshold",
        "essential_threshold",
        "inhibitory_threshold",
        "max_conditions_per_gene",
        "default_search_limit",
        "max_search_limit",
    )

    def __init__(self) -> None:
        self.significance_threshold: float = (
            2.0  # |fitness| > threshold for significant effects
        )
        self.essential_threshold: float = (
            0.5  # Minimum positive fitness for essential genes
        )
        self.inhibitory_threshold: float = (
            -0.5
        )  # Maximum negative fitness for inhibitory genes
        self.max_conditions_per_gene: int = (
            3  # Limit conditions shown per gene in results
        )
        self.default_search_limit: int = 3  # Default limit for search operations
        self.max_search_limit: int = 50  # Maximum allowed search limit


# Global configuration instance
config = FitnessConfig()


# INPUT VALIDATION SECTION
class ValidationResult:
    """Result of input validation for POSIX-style error handling."""

    __slots__ = ("is_valid", "error_message", "suggestions")

    def __init__(
        self,
        is_valid: bool,
        error_message: str = "",
        suggestions: Optional[List[str]] = None,
    ) -> None:
        self.is_valid = is_valid
        self.error_message = error_message
        self.suggestions = suggestions if suggestions is not None else []


def validate_gene_id(gene_id: str) -> ValidationResult:
    """Validate gene ID format before expensive operations.

    POSIX-style: fail fast with clear error messages.
    """
    if not gene_id or not isinstance(gene_id, str):
        return ValidationResult(
            is_valid=False,
            error_message="Gene ID must be a non-empty string",
            suggestions=["Use format like 'Atu0001' or gene symbol like 'rpoA'"],
        )

    if len(gene_id.strip()) == 0:
        return ValidationResult(
            is_valid=False,
            error_message="Gene ID cannot be empty or whitespace",
            suggestions=["Use format like 'Atu0001' or gene symbol like 'rpoA'"],
        )

    return ValidationResult(is_valid=True)


def validate_search_limit(limit: int) -> ValidationResult:
    """Validate search limit parameter for consistency."""
    if not isinstance(limit, int):
        return ValidationResult(
            is_valid=False,
            error_message="Search limit must be an integer",
            suggestions=[f"Use values between 1 and {config.max_search_limit}"],
        )

    if limit < 1:
        return ValidationResult(
            is_valid=False,
            error_message="Search limit must be at least 1",
            suggestions=[f"Use values between 1 and {config.max_search_limit}"],
        )

    if limit > config.max_search_limit:
        return ValidationResult(
            is_valid=False,
            error_message=f"Search limit exceeds maximum of {config.max_search_limit}",
            suggestions=[f"Use values between 1 and {config.max_search_limit}"],
        )

    return ValidationResult(is_valid=True)


def validate_fitness_threshold(
    threshold: float, threshold_type: str = "fitness"
) -> ValidationResult:
    """Validate fitness threshold parameters."""
    if not isinstance(threshold, (int, float)):
        return ValidationResult(
            is_valid=False,
            error_message=f"{threshold_type} threshold must be a number",
            suggestions=["Use numeric values like 0.5, -0.5, or 2.0"],
        )

    return ValidationResult(is_valid=True)


# METADATA REGISTRY SECTION


@dataclass
class GeneMetadata:
    """Metadata for a gene/locus (no fitness data).

    Memory-optimized with __slots__ for better performance.
    """

    __slots__ = ("locus_id", "sys_name", "description")

    locus_id: str
    sys_name: str
    description: str


@dataclass
class ConditionMetadata:
    """Metadata for an experimental condition.

    Memory-optimized with __slots__ for better performance.
    """

    __slots__ = (
        "condition_id",
        "short_desc",
        "long_desc",
        "media",
        "temperature",
        "pH",
        "aerobic",
        "condition_1",
        "concentration_1",
        "units_1",
        "exp_group",
    )

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
    """Metadata for a functional module.

    Memory-optimized with __slots__ for better performance.
    """

    __slots__ = ("module_id", "name", "category", "gene_count", "gene_list")

    module_id: int
    name: str
    category: str
    gene_count: int
    gene_list: List[str]


@dataclass
class FitnessEffect:
    """A pre-filtered significant fitness effect.

    Memory-optimized with __slots__ for better performance.
    """

    __slots__ = ("gene_id", "condition_id", "fitness_value")

    gene_id: str
    condition_id: str
    fitness_value: float


class MetadataRegistry:
    """Centralized registry for all gene, condition, and module metadata.

    Provides efficient, unified access to metadata while keeping fitness data separate.
    Replaces the old ModuleDataLoader and PairsDataLoader with optimized in-memory data structures.
    The FitnessDataLoader is retained for functions requiring complete fitness matrices.
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

        # File tracking
        self.loaded = False
        self._file_mtimes: Dict[str, float] = {}

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
        """Load all metadata from files."""
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


# DATA LOADING SECTION
# Note: ModuleDataLoader and PairsDataLoader have been removed as their functionality
# is now provided by the centralized MetadataRegistry. FitnessDataLoader remains for
# functions that need access to complete fitness matrices rather than just significant effects.


# Legacy ModuleDataLoader class removed - functionality replaced by MetadataRegistry


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

    def _needs_reload(self) -> bool:
        """Check if data needs to be loaded."""
        return not self.loaded

    def load_data(self) -> None:
        """Load the fitness data from the tab-separated file."""
        if not self._needs_reload():
            return

        if not os.path.exists(self.data_file):
            raise FileNotFoundError(f"Fitness data file not found: {self.data_file}")

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

    def search_genes(self, query: str, limit: int = 10) -> List[Dict[str, Any]]:
        """Search for genes by name or description.

        Args:
            query: Search term to match against gene names or descriptions
            limit: Maximum number of results to return

        Returns:
            List of matching gene information
        """
        self.load_data()
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

    def get_conditions(self, condition_filter: Optional[str] = None) -> List[str]:
        """Get list of available growth conditions.

        Args:
            condition_filter: Optional filter to match condition names

        Returns:
            List of condition names
        """
        self.load_data()
        if condition_filter:
            condition_filter = condition_filter.lower()
            return [
                cond for cond in self.conditions if condition_filter in cond.lower()
            ]
        else:
            return self.conditions.copy()

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


# Legacy PairsDataLoader class removed - functionality replaced by MetadataRegistry


# Global data loader instances
fitness_loader = FitnessDataLoader()
# module_loader and pairs_loader removed - functionality replaced by MetadataRegistry


# MCP TOOL SECTION
def get_gene_info(gene_id: str) -> Dict[str, Any]:
    """
    Get basic information about a gene.

    POSIX-style tool: does one thing well - retrieves gene metadata.

    CHAIN WITH:
    - get_gene_fitness() to get fitness data for this gene
    - get_gene_modules() to find functional modules containing this gene
    - search_genes() to find similar genes

    Args:
        gene_id: Gene locus ID (e.g., 'Atu0001') or system name

    Returns:
        Dict with standardized structure:
        {
            "data": {"locusId": str, "sysName": str, "description": str},
            "metadata": {"source": str},
            "suggestions": [str]
        }
    """
    # Input validation
    validation = validate_gene_id(gene_id)
    if not validation.is_valid:
        return {
            "error": validation.error_message,
            "suggestions": validation.suggestions,
        }

    gene_meta = metadata_registry.get_gene(gene_id)
    if gene_meta:
        return {
            "data": {
                "locusId": gene_meta.locus_id,
                "sysName": gene_meta.sys_name,
                "description": gene_meta.description,
            },
            "metadata": {"source": "MetadataRegistry gene lookup"},
            "suggestions": [
                "get_gene_fitness",
                "get_gene_modules",
                "get_fitness_effects_for_gene",
            ],
        }
    else:
        return {
            "error": f"Gene {gene_id} not found",
            "suggestions": [
                "try search_genes() to find similar gene names",
                "verify gene ID format (e.g., 'Atu0001' or gene symbol)",
            ],
        }


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


def search_genes(query: str, limit: int = 3) -> Dict[str, Any]:
    """
    Search for genes by name or description.

    POSIX-style tool: does one thing well - searches gene metadata.

    CHAIN WITH:
    - get_gene_info() to get detailed info for specific genes
    - get_gene_fitness() to analyze fitness data for found genes
    - get_gene_modules() to find functional context

    Args:
        query: Search term to match against gene names or descriptions
        limit: Maximum number of results to return (default: 3)

    Returns:
        Dict with standardized structure:
        {
            "data": {"genes": List[Dict], "total_results": int, "query": str},
            "metadata": {"search_type": str, "limit_applied": int},
            "suggestions": [str]
        }
    """
    # Input validation
    if not query or not isinstance(query, str) or len(query.strip()) == 0:
        return {
            "error": "Search query must be a non-empty string",
            "suggestions": ["Use gene names, symbols, or description keywords"],
        }

    limit_validation = validate_search_limit(limit)
    if not limit_validation.is_valid:
        return {
            "error": limit_validation.error_message,
            "suggestions": limit_validation.suggestions,
        }

    gene_metas = metadata_registry.search_genes(query, limit)

    genes_data = [
        {
            "locusId": gene_meta.locus_id,
            "sysName": gene_meta.sys_name,
            "description": gene_meta.description,
        }
        for gene_meta in gene_metas
    ]

    return {
        "data": {"genes": genes_data, "total_results": len(genes_data), "query": query},
        "metadata": {"search_type": "fuzzy_match", "limit_applied": limit},
        "suggestions": ["get_gene_info", "get_gene_fitness", "get_gene_modules"]
        if genes_data
        else ["try broader search terms", "check spelling of gene names"],
    }


def get_growth_conditions(condition_filter: Optional[str] = None) -> Dict[str, Any]:
    """
    Get list of available growth conditions.

    POSIX-style tool: does one thing well - lists experimental conditions.

    CHAIN WITH:
    - get_condition_details() to get detailed info for specific conditions
    - get_genes_with_fitness_effects() to find genes affected by conditions
    - find_essential_genes() or find_growth_inhibitor_genes() for condition-specific analysis

    Args:
        condition_filter: Optional filter to match condition names (e.g., 'LB', 'glucose', 'metal')

    Returns:
        Dict with standardized structure:
        {
            "data": {"conditions": List[str], "total_conditions": int, "filter_applied": str},
            "metadata": {"source": str},
            "suggestions": [str]
        }
    """
    # Input validation for filter
    if condition_filter is not None and (
        not isinstance(condition_filter, str) or len(condition_filter.strip()) == 0
    ):
        return {
            "error": "Condition filter must be a non-empty string if provided",
            "suggestions": [
                "Use terms like 'pH', 'stress', 'carbon', 'metal', or 'antibiotic'"
            ],
        }

    conditions = metadata_registry.get_all_conditions(condition_filter)

    return {
        "data": {
            "conditions": conditions,
            "total_conditions": len(conditions),
            "filter_applied": condition_filter if condition_filter else "none",
        },
        "metadata": {"source": "MetadataRegistry condition listing"},
        "suggestions": [
            "get_condition_details",
            "get_genes_with_fitness_effects",
            "find_essential_genes",
            "find_growth_inhibitor_genes",
        ]
        if conditions
        else ["try different filter terms", "remove filter to see all conditions"],
    }


def get_condition_details(condition_name: str) -> Dict[str, Any]:
    """
    Get detailed information about a specific growth condition.

    POSIX-style tool: does one thing well - retrieves condition metadata.

    CHAIN WITH:
    - get_genes_with_fitness_effects() to find genes affected by this condition
    - get_growth_conditions() to explore related conditions
    - find_essential_genes() or find_growth_inhibitor_genes() for condition-specific analysis

    Args:
        condition_name: Name of the condition (e.g., 'set10IT004', 'set10IT020')

    Returns:
        Dict with standardized structure:
        {
            "data": {"condition_info": Dict, "experimental_setup": Dict},
            "metadata": {"source": str},
            "suggestions": [str]
        }
    """
    # Input validation
    if (
        not condition_name
        or not isinstance(condition_name, str)
        or len(condition_name.strip()) == 0
    ):
        return {
            "error": "Condition name must be a non-empty string",
            "suggestions": ["Use condition identifiers from get_growth_conditions()"],
        }

    condition_meta = metadata_registry.get_condition(condition_name)
    if not condition_meta:
        return {
            "error": f"Condition {condition_name} not found",
            "suggestions": [
                "try get_growth_conditions() to find valid condition names",
                "check condition name spelling and format",
            ],
        }

    return {
        "data": {
            "condition_info": {
                "condition_name": condition_meta.condition_id,
                "short_description": condition_meta.short_desc,
                "long_description": condition_meta.long_desc,
                "experimental_group": condition_meta.exp_group,
            },
            "experimental_setup": {
                "media": condition_meta.media,
                "temperature": str(condition_meta.temperature) + "°C"
                if condition_meta.temperature
                else "",
                "pH": condition_meta.pH,
                "aerobic": condition_meta.aerobic,
                "treatment": {
                    "compound": condition_meta.condition_1,
                    "concentration": condition_meta.concentration_1,
                    "units": condition_meta.units_1,
                }
                if condition_meta.condition_1
                else None,
            },
        },
        "metadata": {"source": "MetadataRegistry condition lookup"},
        "suggestions": [
            "get_genes_with_fitness_effects",
            "find_essential_genes",
            "find_growth_inhibitor_genes",
            "get_growth_conditions",
        ],
    }


def interpret_fitness_score(fitness_score: float) -> Dict[str, Any]:
    """
    Interpret a fitness score in biological terms.

    POSIX-style tool: does one thing well - interprets fitness values.

    CHAIN WITH:
    - get_fitness_effects_for_gene() to find fitness scores to interpret
    - find_essential_genes() or find_growth_inhibitor_genes() for systematic analysis
    - analyze_gene_fitness() for comprehensive gene analysis

    Args:
        fitness_score: Numerical fitness score from gene knockout experiment

    Returns:
        Dict with standardized structure:
        {
            "data": {"interpretation": str, "effect": str, "magnitude": str, "score": float},
            "metadata": {"interpretation_rules": Dict},
            "suggestions": [str]
        }
    """
    # Input validation
    fitness_validation = validate_fitness_threshold(fitness_score, "fitness score")
    if not fitness_validation.is_valid:
        return {
            "error": fitness_validation.error_message,
            "suggestions": fitness_validation.suggestions,
        }

    interpretation_result = metadata_registry.interpret_fitness_score(fitness_score)

    return {
        "data": {
            "interpretation": interpretation_result["interpretation"],
            "effect": interpretation_result["effect"],
            "magnitude": interpretation_result["magnitude"],
            "score": interpretation_result["score"],
        },
        "metadata": {
            "interpretation_rules": {
                "negative_threshold": -0.1,
                "positive_threshold": 0.1,
                "magnitude_ranges": {
                    "minimal": "<0.2",
                    "moderate": "0.2-0.5",
                    "strong": "0.5-1.0",
                    "very_strong": ">1.0",
                },
            }
        },
        "suggestions": [
            "get_fitness_effects_for_gene",
            "find_essential_genes",
            "find_growth_inhibitor_genes",
        ]
        if interpretation_result["effect"] != "unknown"
        else ["provide a valid numerical fitness score"],
    }


def find_essential_genes(
    condition_filter: Optional[str] = None,
    min_fitness_threshold: float = 0.5,
    limit: int = 5,
) -> Dict[str, Any]:
    """
    Find genes that appear essential (positive fitness scores when knocked out).

    POSIX-style tool: does one thing well - identifies essential genes.

    CHAIN WITH:
    - get_gene_info() to get details for essential genes found
    - interpret_fitness_score() to understand significance levels
    - get_condition_details() to understand experimental conditions
    - analyze_gene_fitness() for detailed analysis of specific genes

    OPTIMIZED: Uses pre-filtered significant fitness effects instead of scanning 4.9M cells.
    Performance: ~40K filtered pairs vs 4.9M cell scan (100x faster).

    Args:
        condition_filter: Optional filter to match specific conditions
        min_fitness_threshold: Minimum fitness score to consider essential (default: 0.5)
        limit: Maximum number of genes to return

    Returns:
        Dict with standardized structure:
        {
            "data": {"essential_genes": List[Dict], "total_genes": int, "analysis_params": Dict},
            "metadata": {"algorithm": str, "performance_info": str},
            "suggestions": [str]
        }
    """
    # Input validation
    if condition_filter is not None and (
        not isinstance(condition_filter, str) or len(condition_filter.strip()) == 0
    ):
        return {
            "error": "Condition filter must be a non-empty string if provided",
            "suggestions": ["Use terms like 'pH', 'stress', 'carbon', 'metal'"],
        }

    threshold_validation = validate_fitness_threshold(
        min_fitness_threshold, "minimum fitness threshold"
    )
    if not threshold_validation.is_valid:
        return {
            "error": threshold_validation.error_message,
            "suggestions": threshold_validation.suggestions,
        }

    limit_validation = validate_search_limit(limit)
    if not limit_validation.is_valid:
        return {
            "error": limit_validation.error_message,
            "suggestions": limit_validation.suggestions,
        }

    metadata_registry.load_data()

    # Use pre-filtered significant effects instead of full matrix scan
    essential_gene_effects: Dict[str, List[FitnessEffect]] = {}

    # Filter significant effects for positive values above threshold
    for effect in metadata_registry.fitness_effects:
        # Apply condition filter if specified
        if (
            condition_filter
            and condition_filter.lower() not in effect.condition_id.lower()
        ):
            continue

        # Filter for essential genes (positive fitness scores)
        if effect.fitness_value >= min_fitness_threshold:
            if effect.gene_id not in essential_gene_effects:
                essential_gene_effects[effect.gene_id] = []
            essential_gene_effects[effect.gene_id].append(effect)

    # Build result structure
    essential_genes = []

    for gene_id, effects in essential_gene_effects.items():
        # Get gene metadata
        gene_meta = metadata_registry.get_gene(gene_id)
        if not gene_meta:
            continue

        # Build condition data with enriched information
        essential_conditions = []
        for effect in effects:
            condition_meta = metadata_registry.get_condition(effect.condition_id)
            essential_conditions.append(
                {
                    "condition": effect.condition_id,
                    "fitness_score": effect.fitness_value,
                    "description": condition_meta.short_desc if condition_meta else "",
                    "interpretation": metadata_registry.interpret_fitness_score(
                        effect.fitness_value
                    ),
                }
            )

        # Sort by fitness score and limit conditions per gene
        essential_conditions.sort(key=itemgetter("fitness_score"), reverse=True)
        essential_conditions = essential_conditions[:3]  # Max 3 conditions per gene

        essential_genes.append(
            {
                "gene": {
                    "locusId": gene_meta.locus_id,
                    "sysName": gene_meta.sys_name,
                    "description": gene_meta.description,
                },
                "essential_in_conditions": essential_conditions,
                "num_essential_conditions": len(essential_conditions),
            }
        )

        if len(essential_genes) >= limit:
            break

    # Sort by number of conditions where gene appears essential
    essential_genes.sort(key=itemgetter("num_essential_conditions"), reverse=True)
    final_genes = essential_genes[:limit]

    return {
        "data": {
            "essential_genes": final_genes,
            "total_genes": len(final_genes),
            "analysis_params": {
                "condition_filter": condition_filter,
                "min_fitness_threshold": min_fitness_threshold,
                "limit": limit,
            },
        },
        "metadata": {
            "algorithm": "pre-filtered significant fitness effects scan",
            "performance_info": "~40K filtered effects vs 4.9M cell scan (100x faster)",
        },
        "suggestions": [
            "get_gene_info",
            "interpret_fitness_score",
            "get_condition_details",
            "analyze_gene_fitness",
        ]
        if final_genes
        else [
            "try lower min_fitness_threshold",
            "try different condition_filter",
            "try find_growth_inhibitor_genes",
        ],
    }


def find_growth_inhibitor_genes(
    condition_filter: Optional[str] = None,
    max_fitness_threshold: float = -0.5,
    limit: int = 5,
) -> Dict[str, Any]:
    """
    Find genes that inhibit growth (negative fitness scores when knocked out).

    POSIX-style tool: does one thing well - identifies growth-inhibiting genes.

    CHAIN WITH:
    - get_gene_info() to get details for inhibitory genes found
    - interpret_fitness_score() to understand significance levels
    - get_condition_details() to understand experimental conditions
    - analyze_gene_fitness() for detailed analysis of specific genes

    OPTIMIZED: Uses pre-filtered significant fitness effects instead of scanning 4.9M cells.
    Performance: ~40K filtered pairs vs 4.9M cell scan (100x faster).

    Args:
        condition_filter: Optional filter to match specific conditions
        max_fitness_threshold: Maximum fitness score to consider inhibitory (default: -0.5)
        limit: Maximum number of genes to return

    Returns:
        Dict with standardized structure:
        {
            "data": {"inhibitory_genes": List[Dict], "total_genes": int, "analysis_params": Dict},
            "metadata": {"algorithm": str, "performance_info": str},
            "suggestions": [str]
        }
    """
    # Input validation
    if condition_filter is not None and (
        not isinstance(condition_filter, str) or len(condition_filter.strip()) == 0
    ):
        return {
            "error": "Condition filter must be a non-empty string if provided",
            "suggestions": ["Use terms like 'pH', 'stress', 'carbon', 'metal'"],
        }

    threshold_validation = validate_fitness_threshold(
        max_fitness_threshold, "maximum fitness threshold"
    )
    if not threshold_validation.is_valid:
        return {
            "error": threshold_validation.error_message,
            "suggestions": threshold_validation.suggestions,
        }

    limit_validation = validate_search_limit(limit)
    if not limit_validation.is_valid:
        return {
            "error": limit_validation.error_message,
            "suggestions": limit_validation.suggestions,
        }

    metadata_registry.load_data()

    # Use pre-filtered significant effects instead of full matrix scan
    inhibitor_gene_effects: Dict[str, List[FitnessEffect]] = {}

    # Filter significant effects for negative values below threshold
    for effect in metadata_registry.fitness_effects:
        # Apply condition filter if specified
        if (
            condition_filter
            and condition_filter.lower() not in effect.condition_id.lower()
        ):
            continue

        # Filter for inhibitory genes (negative fitness scores)
        if effect.fitness_value <= max_fitness_threshold:
            if effect.gene_id not in inhibitor_gene_effects:
                inhibitor_gene_effects[effect.gene_id] = []
            inhibitor_gene_effects[effect.gene_id].append(effect)

    # Build result structure
    inhibitor_genes = []

    for gene_id, effects in inhibitor_gene_effects.items():
        # Get gene metadata
        gene_meta = metadata_registry.get_gene(gene_id)
        if not gene_meta:
            continue

        # Build condition data with enriched information
        inhibitor_conditions = []
        for effect in effects:
            condition_meta = metadata_registry.get_condition(effect.condition_id)
            inhibitor_conditions.append(
                {
                    "condition": effect.condition_id,
                    "fitness_score": effect.fitness_value,
                    "description": condition_meta.short_desc if condition_meta else "",
                    "interpretation": metadata_registry.interpret_fitness_score(
                        effect.fitness_value
                    ),
                }
            )

        # Sort by fitness score and limit conditions per gene
        inhibitor_conditions.sort(key=itemgetter("fitness_score"))
        inhibitor_conditions = inhibitor_conditions[:3]  # Max 3 conditions per gene

        inhibitor_genes.append(
            {
                "gene": {
                    "locusId": gene_meta.locus_id,
                    "sysName": gene_meta.sys_name,
                    "description": gene_meta.description,
                },
                "inhibits_growth_in_conditions": inhibitor_conditions,
                "num_inhibitory_conditions": len(inhibitor_conditions),
            }
        )

        if len(inhibitor_genes) >= limit:
            break

    # Sort by number of conditions where gene inhibits growth
    inhibitor_genes.sort(key=itemgetter("num_inhibitory_conditions"), reverse=True)
    final_genes = inhibitor_genes[:limit]

    return {
        "data": {
            "inhibitory_genes": final_genes,
            "total_genes": len(final_genes),
            "analysis_params": {
                "condition_filter": condition_filter,
                "max_fitness_threshold": max_fitness_threshold,
                "limit": limit,
            },
        },
        "metadata": {
            "algorithm": "pre-filtered significant fitness effects scan",
            "performance_info": "~40K filtered effects vs 4.9M cell scan (100x faster)",
        },
        "suggestions": [
            "get_gene_info",
            "interpret_fitness_score",
            "get_condition_details",
            "analyze_gene_fitness",
        ]
        if final_genes
        else [
            "try higher max_fitness_threshold (less negative)",
            "try different condition_filter",
            "try find_essential_genes",
        ],
    }


def analyze_gene_fitness(
    gene_id: str,
    min_fitness: Optional[float] = None,
    max_fitness: Optional[float] = None,
    limit: Optional[int] = 10,
) -> Dict[str, Any]:
    """
    Analyze fitness effects for a gene knockout mutant across conditions.

    POSIX-style tool: does one thing well - comprehensive gene fitness analysis.

    CHAIN WITH:
    - get_gene_info() to get basic gene information
    - interpret_fitness_score() to understand individual fitness values
    - get_condition_details() to understand experimental conditions
    - find_essential_genes() or find_growth_inhibitor_genes() for focused analysis

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

    POSIX-style tool: does one thing well - retrieves functional modules for a gene.

    CHAIN WITH:
    - get_module_genes() to explore other genes in the same modules
    - search_modules() to find related functional modules
    - get_gene_info() to get basic gene information

    Args:
        gene_id: Gene locus tag (e.g., 'Atu0001')

    Returns:
        Dict with standardized structure:
        {
            "data": {"gene_id": str, "modules": List[Dict], "module_count": int},
            "metadata": {"source": str},
            "suggestions": [str]
        }
    """
    # Input validation
    validation = validate_gene_id(gene_id)
    if not validation.is_valid:
        return {
            "error": validation.error_message,
            "suggestions": validation.suggestions,
        }

    module_metas = metadata_registry.get_gene_modules(gene_id)

    if not module_metas:
        return {
            "error": f"No modules found for gene {gene_id}",
            "suggestions": [
                "try get_gene_info() to verify gene exists",
                "try search_modules() to explore available modules",
            ],
        }

    modules_data = [
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

    return {
        "data": {
            "gene_id": gene_id,
            "modules": modules_data,
            "module_count": len(modules_data),
        },
        "metadata": {"source": "MetadataRegistry module assignments"},
        "suggestions": ["get_module_genes", "search_modules", "get_gene_info"],
    }


def get_module_genes(module_id: int) -> Dict[str, Any]:
    """
    Get all genes in a specific module.

    POSIX-style tool: does one thing well - retrieves genes in a functional module.

    CHAIN WITH:
    - get_gene_info() to get details for genes in the module
    - get_gene_modules() to explore other modules for these genes
    - search_modules() to find related functional modules

    Args:
        module_id: Module ID number

    Returns:
        Dict with standardized structure:
        {
            "data": {"module": Dict, "genes": List[Dict], "gene_count": int},
            "metadata": {"source": str},
            "suggestions": [str]
        }
    """
    # Input validation
    if not isinstance(module_id, int):
        return {
            "error": "Module ID must be an integer",
            "suggestions": [
                "Use module IDs from get_all_modules() or search_modules()"
            ],
        }

    module_meta = metadata_registry.get_module(module_id)
    if not module_meta:
        return {
            "error": f"Module {module_id} not found",
            "suggestions": [
                "try get_all_modules() to find valid module IDs",
                "try search_modules() to search by name or category",
            ],
        }

    gene_metas = metadata_registry.get_module_genes(module_id)

    genes_data = [
        {
            "locus_tag": gene_meta.locus_id,
            "module_id": module_id,
            "gene_weight": 1.0,  # Default weight
            "product": gene_meta.description,
        }
        for gene_meta in gene_metas
    ]

    return {
        "data": {
            "module": {
                "module_id": module_meta.module_id,
                "name": module_meta.name,
                "category": module_meta.category,
                "count": module_meta.gene_count,
            },
            "genes": genes_data,
            "gene_count": len(genes_data),
        },
        "metadata": {"source": "MetadataRegistry module-gene assignments"},
        "suggestions": ["get_gene_info", "get_gene_modules", "search_modules"],
    }


def search_modules(query: str, limit: int = 3) -> Dict[str, Any]:
    """
    Search for modules by name or category.

    POSIX-style tool: does one thing well - searches functional modules.

    CHAIN WITH:
    - get_module_genes() to explore genes in found modules
    - get_gene_modules() to find modules for specific genes
    - get_all_modules() to browse all available modules

    Args:
        query: Search term to match against module names or categories
        limit: Maximum number of results to return (default: 3)

    Returns:
        Dict with standardized structure:
        {
            "data": {"modules": List[Dict], "total_results": int, "query": str},
            "metadata": {"search_type": str, "limit_applied": int},
            "suggestions": [str]
        }
    """
    # Input validation
    if not query or not isinstance(query, str) or len(query.strip()) == 0:
        return {
            "error": "Search query must be a non-empty string",
            "suggestions": ["Use module names, categories, or functional keywords"],
        }

    limit_validation = validate_search_limit(limit)
    if not limit_validation.is_valid:
        return {
            "error": limit_validation.error_message,
            "suggestions": limit_validation.suggestions,
        }

    module_metas = metadata_registry.search_modules(query, limit)

    modules_data = []
    for module_meta in module_metas:
        gene_metas = metadata_registry.get_module_genes(module_meta.module_id)

        genes_data = [
            {
                "locus_tag": gene_meta.locus_id,
                "module_id": module_meta.module_id,
                "gene_weight": 1.0,
                "product": gene_meta.description,
            }
            for gene_meta in gene_metas
        ]

        modules_data.append(
            {
                "module": {
                    "module_id": module_meta.module_id,
                    "name": module_meta.name,
                    "category": module_meta.category,
                    "count": module_meta.gene_count,
                },
                "genes": genes_data,
                "gene_count": len(genes_data),
            }
        )

    return {
        "data": {
            "modules": modules_data,
            "total_results": len(modules_data),
            "query": query,
        },
        "metadata": {"search_type": "fuzzy_match", "limit_applied": limit},
        "suggestions": ["get_module_genes", "get_gene_modules", "get_all_modules"]
        if modules_data
        else [
            "try broader search terms",
            "try get_all_modules() to browse all modules",
        ],
    }


def get_all_modules() -> Dict[str, Any]:
    """
    Get list of all available modules.

    POSIX-style tool: does one thing well - lists all functional modules.

    CHAIN WITH:
    - get_module_genes() to explore genes in specific modules
    - search_modules() to search for modules by keyword
    - get_gene_modules() to find modules for specific genes

    Returns:
        Dict with standardized structure:
        {
            "data": {"modules": List[Dict], "total_modules": int},
            "metadata": {"source": str},
            "suggestions": [str]
        }
    """
    modules_list = metadata_registry.get_all_modules_list()

    return {
        "data": {"modules": modules_list, "total_modules": len(modules_list)},
        "metadata": {"source": "MetadataRegistry complete module listing"},
        "suggestions": ["get_module_genes", "search_modules", "get_gene_modules"]
        if modules_list
        else ["no modules available in current dataset"],
    }


def get_fitness_effects_for_gene(gene_id: str) -> Dict[str, Any]:
    """
    Get all fitness effects for a gene (conditions where |fitness| > threshold).

    POSIX-style tool: does one thing well - retrieves significant fitness effects.

    CHAIN WITH:
    - expand_fitness_network() to explore related genes/conditions
    - interpret_fitness_score() to understand biological meaning
    - get_gene_info() for basic gene metadata

    Args:
        gene_id: Gene locus ID (e.g., 'Atu0001') or gene symbol (e.g., 'rpoA')

    Returns:
        Dict with standardized structure:
        {
            "data": {"gene_id": str, "fitness_effects": List[Dict], "total_effects": int},
            "metadata": {"threshold": float, "interpretation": str},
            "suggestions": ["expand_fitness_network", "interpret_fitness_score"]
        }
    """
    # Input validation
    validation = validate_gene_id(gene_id)
    if not validation.is_valid:
        return {
            "error": validation.error_message,
            "suggestions": validation.suggestions,
        }
    fitness_effects = metadata_registry.get_gene_fitness_effects(gene_id)

    if not fitness_effects:
        return {
            "error": f"No significant fitness effects found for gene {gene_id}",
            "suggestions": [
                "try search_genes() to find similar genes",
                "try analyze_gene_fitness() for complete fitness data",
            ],
        }

    # Convert to standardized format and sort by absolute value of fitness score
    effects_data = [
        {"condition": effect.condition_id, "fitness_value": effect.fitness_value}
        for effect in fitness_effects
    ]
    effects_data.sort(
        key=lambda x: abs(x["fitness_value"])
        if isinstance(x["fitness_value"], (int, float))
        else 0,
        reverse=True,
    )

    return {
        "data": {
            "gene_id": gene_id,
            "fitness_effects": effects_data,
            "total_effects": len(effects_data),
        },
        "metadata": {
            "threshold": config.significance_threshold,
            "interpretation": f"Conditions where gene knockout has significant fitness effects (|value| > {config.significance_threshold})",
        },
        "suggestions": [
            "expand_fitness_network",
            "interpret_fitness_score",
            "get_gene_info",
        ],
    }


def get_genes_with_fitness_effects(condition_id: str) -> Dict[str, Any]:
    """
    Get all genes with significant fitness effects for a specific condition.

    POSIX-style tool: does one thing well - retrieves genes with significant fitness effects.

    CHAIN WITH:
    - expand_fitness_network() to explore related genes/conditions
    - interpret_fitness_score() to understand biological meaning
    - get_condition_details() for experimental condition metadata

    Args:
        condition_id: Condition identifier (e.g., 'set10IT004 D,L-Malic Acid (C)')

    Returns:
        Dict with standardized structure:
        {
            "data": {"condition_id": str, "fitness_effects": List[Dict], "total_genes": int},
            "metadata": {"threshold": float, "interpretation": str},
            "suggestions": [str]
        }
    """
    # Input validation
    if (
        not condition_id
        or not isinstance(condition_id, str)
        or len(condition_id.strip()) == 0
    ):
        return {
            "error": "Condition ID must be a non-empty string",
            "suggestions": [
                "Use valid condition identifiers from get_growth_conditions()"
            ],
        }

    fitness_effects = metadata_registry.get_condition_fitness_effects(condition_id)

    if not fitness_effects:
        return {
            "error": f"No genes with significant fitness effects found for condition {condition_id}",
            "suggestions": [
                "try get_growth_conditions() to find valid condition IDs",
                "try get_condition_details() for condition information",
            ],
        }

    # Convert to standardized format and sort by absolute value of fitness score
    genes_data = [
        {"gene": effect.gene_id, "fitness_value": effect.fitness_value}
        for effect in fitness_effects
    ]
    genes_data.sort(
        key=lambda x: abs(x["fitness_value"])
        if isinstance(x["fitness_value"], (int, float))
        else 0,
        reverse=True,
    )

    return {
        "data": {
            "condition_id": condition_id,
            "fitness_effects": genes_data,
            "total_genes": len(genes_data),
        },
        "metadata": {
            "threshold": config.significance_threshold,
            "interpretation": f"Genes where knockout has significant fitness effects (|value| > {config.significance_threshold}) in this condition",
        },
        "suggestions": [
            "expand_fitness_network",
            "interpret_fitness_score",
            "get_condition_details",
        ],
    }


def expand_fitness_network(gene_id: str, condition_id: str) -> Dict[str, Any]:
    """
    Expand a fitness network from a gene-condition pair to discover related biological interactions.

    POSIX-style tool: does one thing well - network expansion from fitness effects.

    CHAIN WITH:
    - get_fitness_effects_for_gene() to find starting points
    - get_genes_with_fitness_effects() to find related genes
    - interpret_fitness_score() to understand biological significance

    Starting from a specific gene-condition pair, performs two-hop expansion:
    1. Finds all conditions where the gene has significant fitness effects
    2. Finds all genes with significant fitness effects in the condition
    3. Expands to find all conditions for the gene set and all genes for the condition set

    Args:
        gene_id: Gene locus ID (e.g., 'Atu0001') or gene symbol (e.g., 'rpoA')
        condition_id: Condition identifier (e.g., 'set10IT004 D,L-Malic Acid (C)')

    Returns:
        Dict with standardized network structure and expansion metrics
    """
    # Input validation
    gene_validation = validate_gene_id(gene_id)
    if not gene_validation.is_valid:
        return {
            "error": gene_validation.error_message,
            "suggestions": gene_validation.suggestions,
        }

    if (
        not condition_id
        or not isinstance(condition_id, str)
        or len(condition_id.strip()) == 0
    ):
        return {
            "error": "Condition ID must be a non-empty string",
            "suggestions": [
                "Use valid condition identifiers from get_growth_conditions()"
            ],
        }

    metadata_registry.load_data()

    # Check if the gene-condition pair exists (use index for quick lookup)
    conditions_for_gene = set(metadata_registry.gene_to_conditions.get(gene_id, []))

    if condition_id not in conditions_for_gene:
        return {
            "error": f"No significant fitness effect found for gene {gene_id} in condition {condition_id}",
            "suggestions": [
                "try get_fitness_effects_for_gene() to find valid gene-condition pairs",
                "try get_genes_with_fitness_effects() to find genes for this condition",
            ],
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
mcp.tool(get_fitness_effects_for_gene)
mcp.tool(get_genes_with_fitness_effects)
mcp.tool(expand_fitness_network)


def main() -> None:
    """Main entry point for the application."""
    mcp.run()


if __name__ == "__main__":
    main()
