################################################################################
# fitness_mcp/main.py
# FastMCP wrapper for Agrobacterium mutant fitness data analysis
#
# This MCP analyzes fitness data from barcoded Agrobacterium mutants grown in
# mixed cultures across different conditions. Each row represents a gene knockout
# mutant, and fitness scores indicate:
# - NEGATIVE values: Gene knockout IMPROVES fitness (gene normally inhibits growth)
# - POSITIVE values: Gene knockout REDUCES fitness (gene is beneficial/essential)
# - Values near 0: Gene knockout has minimal effect on fitness
################################################################################
import csv
import os
from typing import Any, Dict, List, Optional
from fastmcp import FastMCP


# DATA LOADING SECTION
import threading
from functools import lru_cache


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
            - Negative values: Gene knockout improves fitness (gene normally inhibits growth)
            - Positive values: Gene knockout reduces fitness (gene is beneficial/essential)
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
            effect = "gene_inhibits_growth"
            interpretation = f"Gene knockout improves fitness ({magnitude} effect). This gene normally inhibits growth in this condition."
        elif fitness_score > 0.1:
            effect = "gene_benefits_growth"
            interpretation = f"Gene knockout reduces fitness ({magnitude} effect). This gene is beneficial/essential for growth in this condition."
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


# Global data loader instances
fitness_loader = FitnessDataLoader()
module_loader = ModuleDataLoader()


# MCP TOOL SECTION
def get_gene_info(gene_id: str) -> Dict[str, Any]:
    """
    Get basic information about a gene.

    Args:
        gene_id: Gene locus ID (e.g., 'Atu0001') or system name

    Returns:
        Dict containing gene information including locusId, sysName, and description
    """
    gene_info = fitness_loader.get_gene_info(gene_id)
    if gene_info:
        return {
            "locusId": gene_info["locusId"],
            "sysName": gene_info["sysName"],
            "description": gene_info["description"],
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
    return fitness_loader.get_gene_fitness(gene_id, condition_filter)


def search_genes(query: str, limit: int = 10) -> List[Dict[str, Any]]:
    """
    Search for genes by name or description.

    Args:
        query: Search term to match against gene names or descriptions
        limit: Maximum number of results to return (default: 10)

    Returns:
        List of dictionaries containing matching gene information
    """
    return fitness_loader.search_genes(query, limit)


def get_growth_conditions(condition_filter: Optional[str] = None) -> List[str]:
    """
    Get list of available growth conditions.

    Args:
        condition_filter: Optional filter to match condition names (e.g., 'LB', 'glucose', 'metal')

    Returns:
        List of condition names
    """
    return fitness_loader.get_conditions(condition_filter)


def get_condition_details(condition_name: str) -> Dict[str, Any]:
    """
    Get detailed information about a specific growth condition.

    Args:
        condition_name: Name of the condition (e.g., 'set10IT004', 'set10IT020')

    Returns:
        Dict with detailed condition information including description, media, temperature, pH, and treatment details
    """
    return fitness_loader.get_condition_details(condition_name)


def interpret_fitness_score(fitness_score: float) -> Dict[str, Any]:
    """
    Interpret a fitness score in biological terms.

    Args:
        fitness_score: Numerical fitness score from gene knockout experiment

    Returns:
        Dict with biological interpretation of the fitness effect:
        - Negative scores: Gene knockout improves fitness (gene normally inhibits growth)
        - Positive scores: Gene knockout reduces fitness (gene is beneficial/essential)
        - Near zero: Gene knockout has minimal effect
    """
    return fitness_loader.interpret_fitness_score(fitness_score)


def find_essential_genes(
    condition_filter: Optional[str] = None,
    min_fitness_threshold: float = 0.5,
    limit: int = 20,
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
    limit: int = 20,
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
) -> Dict[str, Any]:
    """
    Analyze fitness effects for a gene knockout mutant across conditions.

    Args:
        gene_id: Gene locus ID or system name
        min_fitness: Minimum fitness value to include (use negative values to find where gene inhibits growth)
        max_fitness: Maximum fitness value to include (use positive values to find where gene benefits growth)

    Returns:
        Dict with gene info and categorized fitness data. Categories indicate:
        - conditions_where_gene_inhibits_growth: Negative fitness scores
        - conditions_where_gene_benefits_growth: Positive fitness scores
        - neutral_conditions: Fitness scores near zero
    """
    fitness_data = fitness_loader.get_gene_fitness(gene_id)

    if "error" in fitness_data:
        return fitness_data

    # Categorize conditions by fitness effect (corrected interpretation)
    gene_inhibits_growth = []  # Negative fitness: gene knockout improves fitness
    gene_benefits_growth = []  # Positive fitness: gene knockout reduces fitness
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
        ):  # Gene knockout improves fitness (gene normally inhibits growth)
            gene_inhibits_growth.append(item)
        elif (
            fitness_val > 0.5
        ):  # Gene knockout reduces fitness (gene is beneficial/essential)
            gene_benefits_growth.append(item)
        else:  # Neutral effect
            neutral.append(item)

    return {
        "gene": fitness_data["gene"],
        "analysis": {
            "conditions_where_gene_inhibits_growth": gene_inhibits_growth,
            "conditions_where_gene_benefits_growth": gene_benefits_growth,
            "neutral_conditions": neutral,
            "summary": {
                "total_conditions_tested": len(
                    [
                        x
                        for x in fitness_data["fitness_data"]
                        if x["fitness"] is not None
                    ]
                ),
                "inhibitory_count": len(gene_inhibits_growth),
                "beneficial_count": len(gene_benefits_growth),
                "neutral_count": len(neutral),
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
    modules = module_loader.get_modules_for_gene(gene_id)

    if not modules:
        return {"error": f"No modules found for gene {gene_id}"}

    return {"gene_id": gene_id, "modules": modules, "module_count": len(modules)}


def get_module_genes(module_id: int) -> Dict[str, Any]:
    """
    Get all genes in a specific module.

    Args:
        module_id: Module ID number

    Returns:
        Dict containing module info and all genes in the module
    """
    return module_loader.get_genes_in_module(module_id)


def search_modules(query: str, limit: int = 10) -> List[Dict[str, Any]]:
    """
    Search for modules by name or category.

    Args:
        query: Search term to match against module names or categories
        limit: Maximum number of results to return (default: 10)

    Returns:
        List of matching modules with their genes
    """
    return module_loader.search_modules_by_name(query, limit)


def get_all_modules() -> List[Dict[str, Any]]:
    """
    Get list of all available modules.

    Returns:
        List of all modules with basic information
    """
    return module_loader.get_all_modules()


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


def main() -> None:
    """Main entry point for the application."""
    mcp.run()


if __name__ == "__main__":
    main()
