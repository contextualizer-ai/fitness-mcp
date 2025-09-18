"""Tests for data loading classes in fitness_mcp.main module."""

import csv
import os
import tempfile
import threading
import pytest


from src.fitness_mcp.main import FitnessDataLoader
# Note: ModuleDataLoader was removed - functionality replaced by MetadataRegistry


class TestFitnessDataLoader:
    """Test cases for FitnessDataLoader class."""

    def setup_method(self):
        """Set up test fixtures before each test method."""
        self.temp_dir = tempfile.mkdtemp()
        self.data_file = os.path.join(self.temp_dir, "test_data.tab")
        self.exp_desc_file = os.path.join(self.temp_dir, "test_exp.txt")

        # Sample data for testing
        self.sample_data = [
            ["locusId", "sysName", "desc", "condition1", "condition2", "condition3"],
            ["Atu0001", "geneA", "Gene A description", "1.5", "-2.1", "0.1"],
            ["Atu0002", "geneB", "Gene B description", "-0.8", "3.2", "-1.0"],
            ["Atu0003", "geneC", "Gene C description", "0.2", "0.5", "-2.5"],
        ]

        # Sample experimental description data
        self.sample_exp_data = [
            [
                "expName",
                "expDesc",
                "expDescLong",
                "media",
                "temperature",
                "pH",
                "aerobic",
                "condition_1",
                "concentration_1",
                "units_1",
                "expGroup",
            ],
            [
                "condition1",
                "Stress test",
                "Stress condition long desc",
                "LB",
                "30",
                "7.0",
                "TRUE",
                "heat",
                "42",
                "C",
                "stress",
            ],
            [
                "condition2",
                "Carbon test",
                "Carbon source test",
                "M9",
                "25",
                "7.2",
                "TRUE",
                "glucose",
                "10",
                "mM",
                "carbon",
            ],
        ]

    def teardown_method(self):
        """Clean up after each test method."""
        import shutil

        shutil.rmtree(self.temp_dir)

    def create_test_files(self):
        """Create test data files."""
        # Create main data file
        with open(self.data_file, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(self.sample_data)

        # Create experimental description file
        with open(self.exp_desc_file, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(self.sample_exp_data)

    def test_init(self):
        """Test FitnessDataLoader initialization."""
        loader = FitnessDataLoader(self.data_file, self.exp_desc_file)

        assert loader.data_file == self.data_file
        assert loader.exp_desc_file == self.exp_desc_file
        assert not loader.loaded
        assert len(loader.genes) == 0
        assert len(loader.conditions) == 0

    def test_needs_reload_no_files(self):
        """Test _needs_reload when files don't exist."""
        loader = FitnessDataLoader(self.data_file, self.exp_desc_file)
        assert loader._needs_reload()  # Should return True if data not loaded yet

    def test_needs_reload_with_files(self):
        """Test _needs_reload with existing files."""
        self.create_test_files()
        loader = FitnessDataLoader(self.data_file, self.exp_desc_file)

        # Should need reload initially
        assert loader._needs_reload()

        # Load data
        loader.load_data()

        # Should not need reload immediately after loading
        assert not loader._needs_reload()

    def test_load_data_missing_files(self):
        """Test load_data when data file is missing."""
        loader = FitnessDataLoader(self.data_file, self.exp_desc_file)

        with pytest.raises(FileNotFoundError):
            loader.load_data()

    def test_load_data_success(self):
        """Test successful data loading."""
        self.create_test_files()
        loader = FitnessDataLoader(self.data_file, self.exp_desc_file)

        loader.load_data()

        assert loader.loaded
        assert len(loader.genes) >= 3  # Should have loaded 3 genes
        assert len(loader.conditions) == 3  # Should have 3 conditions
        assert "Atu0001" in loader.genes
        assert "geneA" in loader.genes  # Should be indexed by sys_name too

    def test_load_condition_details(self):
        """Test loading of condition details."""
        self.create_test_files()
        loader = FitnessDataLoader(self.data_file, self.exp_desc_file)

        loader.load_data()

        assert len(loader.condition_details) >= 2
        assert "condition1" in loader.condition_details
        details = loader.condition_details["condition1"]
        assert details["short_desc"] == "Stress test"
        assert details["media"] == "LB"

    def test_get_gene_info(self):
        """Test get_gene_info method."""
        self.create_test_files()
        loader = FitnessDataLoader(self.data_file, self.exp_desc_file)

        # Test with non-existent gene
        result = loader.get_gene_info("nonexistent")
        assert result is None

        # Test with existing gene
        loader.load_data()
        result = loader.get_gene_info("Atu0001")
        assert result is not None
        assert result["locusId"] == "Atu0001"
        assert result["sysName"] == "geneA"

    def test_get_gene_fitness(self):
        """Test get_gene_fitness method."""
        self.create_test_files()
        loader = FitnessDataLoader(self.data_file, self.exp_desc_file)

        # Test with non-existent gene
        result = loader.get_gene_fitness("nonexistent")
        assert "error" in result

        # Test with existing gene
        loader.load_data()
        result = loader.get_gene_fitness("Atu0001")
        assert "gene" in result
        assert "fitness_data" in result
        assert result["gene"]["locusId"] == "Atu0001"
        assert len(result["fitness_data"]) == 3

    def test_get_gene_fitness_with_filter(self):
        """Test get_gene_fitness with condition filter."""
        self.create_test_files()
        loader = FitnessDataLoader(self.data_file, self.exp_desc_file)

        loader.load_data()
        result = loader.get_gene_fitness("Atu0001", condition_filter="condition1")
        assert "fitness_data" in result
        # Should only return conditions matching the filter
        matching_conditions = [
            item for item in result["fitness_data"] if "condition1" in item["condition"]
        ]
        assert len(matching_conditions) >= 1

    def test_search_genes(self):
        """Test search_genes method."""
        self.create_test_files()
        loader = FitnessDataLoader(self.data_file, self.exp_desc_file)

        loader.load_data()

        # Search by description
        results = loader.search_genes("Gene A")
        assert len(results) >= 1
        assert any(gene["locusId"] == "Atu0001" for gene in results)

        # Search by sys_name
        results = loader.search_genes("geneB")
        assert len(results) >= 1
        assert any(gene["locusId"] == "Atu0002" for gene in results)

    def test_get_conditions(self):
        """Test get_conditions method."""
        self.create_test_files()
        loader = FitnessDataLoader(self.data_file, self.exp_desc_file)

        loader.load_data()

        # Get all conditions
        conditions = loader.get_conditions()
        assert len(conditions) == 3

        # Get with filter
        filtered = loader.get_conditions("condition1")
        assert len(filtered) >= 1
        assert "condition1" in filtered

    def test_interpret_fitness_score(self):
        """Test interpret_fitness_score method."""
        loader = FitnessDataLoader()

        # Test strongly negative score (essential gene)
        result = loader.interpret_fitness_score(-2.0)
        assert result["effect"] == "gene_benefits_growth"
        assert "beneficial/essential" in result["interpretation"]

        # Test strongly positive score (inhibitory gene)
        result = loader.interpret_fitness_score(2.0)
        assert result["effect"] == "gene_inhibits_growth"
        assert "inhibits growth" in result["interpretation"]

        # Test neutral score
        result = loader.interpret_fitness_score(0.05)
        assert result["effect"] == "neutral"

        # Test None value
        result = loader.interpret_fitness_score(None)
        assert result["effect"] == "unknown"

    def test_get_condition_details(self):
        """Test get_condition_details method."""
        self.create_test_files()
        loader = FitnessDataLoader(self.data_file, self.exp_desc_file)

        # Test with non-existent condition
        result = loader.get_condition_details("nonexistent")
        assert "error" in result

        # Test with existing condition
        loader.load_data()
        result = loader.get_condition_details("condition1")
        assert "condition_name" in result
        assert result["condition_name"] == "condition1"
        assert "short_description" in result
        assert result["short_description"] == "Stress test"

    def test_thread_safety(self):
        """Test thread safety of load_data method."""
        self.create_test_files()
        loader = FitnessDataLoader(self.data_file, self.exp_desc_file)

        results = []
        exceptions = []

        def load_data_thread():
            try:
                loader.load_data()
                results.append(len(loader.genes))
            except Exception as e:
                exceptions.append(e)

        # Start multiple threads
        threads = []
        for _ in range(5):
            thread = threading.Thread(target=load_data_thread)
            threads.append(thread)
            thread.start()

        # Wait for all threads to complete
        for thread in threads:
            thread.join()

        # Should have no exceptions and consistent results
        assert len(exceptions) == 0
        assert all(
            result == results[0] for result in results
        )  # All results should be consistent
        assert results[0] > 0  # Should have loaded some genes


# TestModuleDataLoader removed - functionality replaced by MetadataRegistry in Issue #22
# Equivalent tests are in tests/test_metadata_registry.py
