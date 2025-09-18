"""Tests for the MetadataRegistry centralized metadata management system."""

import tempfile
import pytest
from pathlib import Path

from src.fitness_mcp.main import (
    MetadataRegistry,
    GeneMetadata,
    ConditionMetadata,
    ModuleMetadata,
    FitnessEffect,
)


class TestMetadataRegistry:
    """Test the centralized MetadataRegistry system."""

    @pytest.fixture
    def temp_data_dir(self):
        """Create temporary directory with mock data files."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Create mock fitness data file
            fitness_file = temp_path / "fit_t.tab"
            fitness_content = """locusId	sysName	desc	condition1	condition2	condition3
Atu0001	geneA	Gene A description	1.5	-2.1	0.1
Atu0002	geneB	Gene B description	-0.8	3.2	-1.0
Atu0003	geneC	Gene C description	0.2	0.5	-2.5"""
            fitness_file.write_text(fitness_content)

            # Create mock experimental conditions file
            exp_file = temp_path / "exp_organism_Agro.txt"
            exp_content = """expName	expDesc	expDescLong	media	temperature	pH	aerobic	condition_1	concentration_1	units_1	expGroup
condition1	Stress test	Stress condition long desc	LB	30	7.0	TRUE	heat	42	C	stress
condition2	Carbon test	Carbon source test	M9	25	7.2	TRUE	glucose	10	mM	carbon
condition3	Metal test	Metal stress test	LB	30	7.0	TRUE	copper	5	uM	metal"""
            exp_file.write_text(exp_content)

            # Create mock modules file
            modules_file = temp_path / "RbTnSeq_modules_t1e-7.csv"
            modules_content = """module,locus_tag,gene_weight,product,Description,Preferred_name,Protein names,Gene Ontology (GO)
1,Atu0001,0.8,Product A,Description A,Preferred A,Protein A,GO:123
1,Atu0002,0.6,Product B,Description B,Preferred B,Protein B,GO:456
2,Atu0003,0.9,Product C,Description C,Preferred C,Protein C,GO:789"""
            modules_file.write_text(modules_content)

            # Create mock module metadata file
            module_meta_file = temp_path / "module_meta.tsv"
            module_meta_content = """module	name	category	count
1	Module Alpha	Transport	2
2	Module Beta	Metabolism	1"""
            module_meta_file.write_text(module_meta_content)

            # Create mock pairs file
            pairs_file = temp_path / "fit_t_pairs_threshold_2_long.tab"
            pairs_content = """gene_id	condition_id	value
Atu0002	condition2	3.2
Atu0003	condition3	-2.5"""
            pairs_file.write_text(pairs_content)

            yield {
                "temp_dir": temp_path,
                "fitness_file": str(fitness_file),
                "exp_file": str(exp_file),
                "modules_file": str(modules_file),
                "module_meta_file": str(module_meta_file),
                "pairs_file": str(pairs_file),
            }

    @pytest.fixture
    def registry(self, temp_data_dir):
        """Create MetadataRegistry with temp data files."""
        return MetadataRegistry(
            fitness_file=temp_data_dir["fitness_file"],
            exp_desc_file=temp_data_dir["exp_file"],
            modules_file=temp_data_dir["modules_file"],
            module_meta_file=temp_data_dir["module_meta_file"],
            pairs_file=temp_data_dir["pairs_file"],
        )

    def test_initialization(self, registry):
        """Test registry initialization."""
        assert not registry.loaded
        assert len(registry.genes) == 0
        assert len(registry.conditions) == 0
        assert len(registry.modules) == 0
        assert len(registry.fitness_effects) == 0

    def test_load_gene_metadata(self, registry):
        """Test loading gene metadata from fitness file."""
        registry.load_data()

        assert registry.loaded
        assert len(registry.genes) >= 3  # Should have at least 3 unique genes

        # Check that genes are indexed by both locus_id and sys_name
        assert "Atu0001" in registry.genes
        assert "geneA" in registry.genes
        assert registry.genes["Atu0001"] == registry.genes["geneA"]

        # Verify gene metadata structure
        gene_meta = registry.genes["Atu0001"]
        assert isinstance(gene_meta, GeneMetadata)
        assert gene_meta.locus_id == "Atu0001"
        assert gene_meta.sys_name == "geneA"
        assert gene_meta.description == "Gene A description"

    def test_load_condition_metadata(self, registry):
        """Test loading condition metadata."""
        registry.load_data()

        assert len(registry.conditions) == 3
        assert "condition1" in registry.conditions

        condition_meta = registry.conditions["condition1"]
        assert isinstance(condition_meta, ConditionMetadata)
        assert condition_meta.condition_id == "condition1"
        assert condition_meta.short_desc == "Stress test"
        assert condition_meta.media == "LB"
        assert condition_meta.temperature == "30"

    def test_load_module_metadata(self, registry):
        """Test loading module metadata and gene-module assignments."""
        registry.load_data()

        assert len(registry.modules) == 2
        assert 1 in registry.modules
        assert 2 in registry.modules

        module_meta = registry.modules[1]
        assert isinstance(module_meta, ModuleMetadata)
        assert module_meta.module_id == 1
        assert module_meta.name == "Module Alpha"
        assert module_meta.category == "Transport"
        assert module_meta.gene_count == 2
        assert "Atu0001" in module_meta.gene_list
        assert "Atu0002" in module_meta.gene_list

    def test_load_fitness_effects(self, registry):
        """Test loading significant fitness effects."""
        registry.load_data()

        assert len(registry.fitness_effects) == 2

        # Check fitness effects structure
        effect = registry.fitness_effects[0]
        assert isinstance(effect, FitnessEffect)
        assert hasattr(effect, "gene_id")
        assert hasattr(effect, "condition_id")
        assert hasattr(effect, "fitness_value")

    def test_build_indexes(self, registry):
        """Test that lookup indexes are built correctly."""
        registry.load_data()

        # Test gene-to-modules index
        assert "Atu0001" in registry.gene_to_modules
        assert 1 in registry.gene_to_modules["Atu0001"]

        # Test module-to-genes index
        assert 1 in registry.module_to_genes
        assert "Atu0001" in registry.module_to_genes[1]
        assert "Atu0002" in registry.module_to_genes[1]

        # Test gene-to-conditions index (from fitness effects)
        assert "Atu0002" in registry.gene_to_conditions
        assert "condition2" in registry.gene_to_conditions["Atu0002"]

    def test_get_gene(self, registry):
        """Test gene lookup by ID."""
        registry.load_data()

        # Test lookup by locus_id
        gene_meta = registry.get_gene("Atu0001")
        assert gene_meta is not None
        assert gene_meta.locus_id == "Atu0001"

        # Test lookup by sys_name
        gene_meta = registry.get_gene("geneA")
        assert gene_meta is not None
        assert gene_meta.locus_id == "Atu0001"

        # Test non-existent gene
        assert registry.get_gene("NONEXISTENT") is None

    def test_get_condition(self, registry):
        """Test condition lookup by ID."""
        registry.load_data()

        condition_meta = registry.get_condition("condition1")
        assert condition_meta is not None
        assert condition_meta.condition_id == "condition1"

        assert registry.get_condition("NONEXISTENT") is None

    def test_get_module(self, registry):
        """Test module lookup by ID."""
        registry.load_data()

        module_meta = registry.get_module(1)
        assert module_meta is not None
        assert module_meta.module_id == 1

        assert registry.get_module(999) is None

    def test_search_genes(self, registry):
        """Test gene search functionality."""
        registry.load_data()

        # Search by description
        results = registry.search_genes("Gene A", limit=5)
        assert len(results) == 1
        assert results[0].locus_id == "Atu0001"

        # Search by sys_name
        results = registry.search_genes("geneB", limit=5)
        assert len(results) == 1
        assert results[0].locus_id == "Atu0002"

    def test_search_modules(self, registry):
        """Test module search functionality."""
        registry.load_data()

        # Search by name
        results = registry.search_modules("Alpha", limit=5)
        assert len(results) == 1
        assert results[0].module_id == 1

        # Search by category
        results = registry.search_modules("Transport", limit=5)
        assert len(results) == 1
        assert results[0].module_id == 1

    def test_get_gene_modules(self, registry):
        """Test getting modules for a gene."""
        registry.load_data()

        modules = registry.get_gene_modules("Atu0001")
        assert len(modules) == 1
        assert modules[0].module_id == 1

    def test_get_module_genes(self, registry):
        """Test getting genes for a module."""
        registry.load_data()

        genes = registry.get_module_genes(1)
        assert len(genes) == 2
        gene_ids = [gene.locus_id for gene in genes]
        assert "Atu0001" in gene_ids
        assert "Atu0002" in gene_ids

    def test_get_gene_fitness_effects(self, registry):
        """Test getting fitness effects for a gene."""
        registry.load_data()

        effects = registry.get_gene_fitness_effects("Atu0002")
        assert len(effects) == 1
        assert effects[0].condition_id == "condition2"
        assert effects[0].fitness_value == 3.2

    def test_get_condition_fitness_effects(self, registry):
        """Test getting fitness effects for a condition."""
        registry.load_data()

        effects = registry.get_condition_fitness_effects("condition2")
        assert len(effects) == 1
        assert effects[0].gene_id == "Atu0002"
        assert effects[0].fitness_value == 3.2

    def test_interpret_fitness_score(self, registry):
        """Test fitness score interpretation."""
        # Test strongly negative score (essential gene)
        result = registry.interpret_fitness_score(-2.0)
        assert result["effect"] == "gene_benefits_growth"
        assert "beneficial/essential" in result["interpretation"]

        # Test strongly positive score (inhibitory gene)
        result = registry.interpret_fitness_score(2.0)
        assert result["effect"] == "gene_inhibits_growth"
        assert "inhibits growth" in result["interpretation"]

        # Test neutral score
        result = registry.interpret_fitness_score(0.05)
        assert result["effect"] == "neutral"

        # Test None value
        result = registry.interpret_fitness_score(None)
        assert result["effect"] == "unknown"

    def test_get_all_conditions(self, registry):
        """Test getting all conditions with optional filtering."""
        registry.load_data()

        # Get all conditions
        all_conditions = registry.get_all_conditions()
        assert len(all_conditions) == 3
        assert "condition1" in all_conditions

        # Test filtering (matches condition ID)
        filtered = registry.get_all_conditions("condition1")
        assert len(filtered) == 1
        assert filtered[0] == "condition1"

    def test_get_all_modules_list(self, registry):
        """Test getting all modules as list."""
        registry.load_data()

        modules_list = registry.get_all_modules_list()
        assert len(modules_list) == 2
        assert all("module_id" in module for module in modules_list)
        assert all("name" in module for module in modules_list)

    def test_file_change_detection(self, registry, temp_data_dir):
        """Test that registry detects file changes and reloads."""
        # Load data initially
        registry.load_data()
        initial_genes_count = len(registry.genes)

        # Modify fitness file
        fitness_file = Path(temp_data_dir["fitness_file"])
        content = fitness_file.read_text()
        content += "\nAtu0004\tgeneD\tGene D description\t0.5\t-1.0\t2.0\n"
        fitness_file.write_text(content)

        # Force reload by clearing loaded flag and checking needs_reload
        registry.loaded = False
        assert registry._needs_reload()

        # Load data again
        registry.load_data()
        new_genes_count = len(registry.genes)

        # Should have more genes now
        assert new_genes_count > initial_genes_count
        assert "Atu0004" in registry.genes

    def test_thread_safety(self, registry):
        """Test thread safety of load_data method."""
        import threading

        load_results = []

        def load_data_thread():
            try:
                registry.load_data()
                load_results.append(True)
            except Exception as e:
                load_results.append(f"Error: {e}")

        # Start multiple threads loading data simultaneously
        threads = []
        for _ in range(5):
            thread = threading.Thread(target=load_data_thread)
            threads.append(thread)
            thread.start()

        # Wait for all threads to complete
        for thread in threads:
            thread.join()

        # All loads should succeed
        assert len(load_results) == 5
        assert all(result is True for result in load_results)
        assert registry.loaded

    def test_missing_files_handling(self):
        """Test handling of missing data files."""
        # Create registry with non-existent files
        registry = MetadataRegistry(
            fitness_file="nonexistent_fitness.tab",
            exp_desc_file="nonexistent_exp.txt",
            modules_file="nonexistent_modules.csv",
            module_meta_file="nonexistent_meta.tsv",
            pairs_file="nonexistent_pairs.tab",
        )

        # Should not crash, should handle missing files gracefully
        registry.load_data()

        # Data structures should be empty but registry should still work
        assert len(registry.genes) == 0
        assert len(registry.conditions) == 0
        assert len(registry.modules) == 0
        assert len(registry.fitness_effects) == 0

        # Methods should return appropriate empty results
        assert registry.get_gene("any_gene") is None
        assert registry.search_genes("any_query") == []
        assert registry.get_all_conditions() == []

    def test_pairs_file_generation(self, registry, temp_data_dir):
        """Test automatic generation of pairs file when missing."""
        # Remove pairs file
        pairs_file = Path(temp_data_dir["pairs_file"])
        pairs_file.unlink()

        # Should generate pairs file during load
        registry.load_data()

        # Pairs file should be recreated
        assert pairs_file.exists()

        # Should have fitness effects loaded
        assert len(registry.fitness_effects) > 0
