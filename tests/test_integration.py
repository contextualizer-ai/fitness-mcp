"""Integration tests for fitness_mcp using real data files."""

import os
import pytest
from unittest.mock import patch

from src.fitness_mcp.main import (
    fitness_loader,
    module_loader,
    get_gene_info,
    get_gene_fitness,
    search_genes,
    find_essential_genes,
    analyze_gene_fitness,
    get_gene_modules,
    search_modules,
)


class TestIntegrationWithRealData:
    """Integration tests using actual data files if available."""

    @pytest.fixture(autouse=True)
    def setup_data_paths(self):
        """Check if real data files exist and set up paths."""
        self.data_dir = "data"
        self.fitness_file = os.path.join(self.data_dir, "fit_t.tab")
        self.exp_file = os.path.join(self.data_dir, "exp_organism_Agro.txt")
        self.modules_file = os.path.join(self.data_dir, "RbTnSeq_modules_t1e-7.csv")
        self.meta_file = os.path.join(self.data_dir, "module_meta.tsv")

        self.has_fitness_data = os.path.exists(self.fitness_file) and os.path.exists(
            self.exp_file
        )
        self.has_module_data = os.path.exists(self.modules_file) and os.path.exists(
            self.meta_file
        )

    def test_fitness_loader_real_data(self):
        """Test fitness loader with real data files."""
        if not self.has_fitness_data:
            pytest.skip("Real fitness data not available")
        # Temporarily override file paths
        original_data_file = fitness_loader.data_file
        original_exp_file = fitness_loader.exp_desc_file

        try:
            fitness_loader.data_file = self.fitness_file
            fitness_loader.exp_desc_file = self.exp_file
            fitness_loader.loaded = False  # Force reload

            fitness_loader.load_data()

            assert fitness_loader.loaded
            assert len(fitness_loader.genes) > 0
            assert len(fitness_loader.conditions) > 0

            # Test that we can get info for a real gene
            first_gene = next(iter(fitness_loader.genes.keys()))
            gene_info = fitness_loader.get_gene_info(first_gene)
            assert gene_info is not None
            assert "locusId" in gene_info

        finally:
            # Restore original paths
            fitness_loader.data_file = original_data_file
            fitness_loader.exp_desc_file = original_exp_file

    def test_module_loader_real_data(self):
        """Test module loader with real data files."""
        if not self.has_module_data:
            pytest.skip("Real module data not available")
        # Temporarily override file paths
        original_modules_file = module_loader.modules_file
        original_meta_file = module_loader.meta_file

        try:
            module_loader.modules_file = self.modules_file
            module_loader.meta_file = self.meta_file
            module_loader.loaded = False  # Force reload

            module_loader.load_data()

            assert module_loader.loaded
            assert len(module_loader.gene_to_modules) > 0
            assert len(module_loader.module_to_genes) > 0
            assert len(module_loader.module_meta) > 0

        finally:
            # Restore original paths
            module_loader.modules_file = original_modules_file
            module_loader.meta_file = original_meta_file

    def test_search_genes_integration(self):
        """Test search_genes integration with mocked data."""
        mock_genes = {
            "Atu0001": {
                "locusId": "Atu0001",
                "sysName": "rpoA",
                "description": "RNA polymerase alpha subunit",
                "fitness_values": [0.5, -0.3, 0.1],
            },
            "rpoA": {  # Same gene indexed by sysName
                "locusId": "Atu0001",
                "sysName": "rpoA",
                "description": "RNA polymerase alpha subunit",
                "fitness_values": [0.5, -0.3, 0.1],
            },
            "Atu0002": {
                "locusId": "Atu0002",
                "sysName": "dnaA",
                "description": "Chromosomal replication initiator protein DnaA",
                "fitness_values": [1.2, 0.8, 0.9],
            },
        }

        with (
            patch.object(fitness_loader, "genes", mock_genes),
            patch.object(fitness_loader, "loaded", True),
        ):
            # Search by gene function
            results = search_genes("polymerase", limit=5)
            assert len(results) >= 1
            assert any(
                "polymerase" in result["description"].lower() for result in results
            )

            # Search by gene name
            results = search_genes("rpoA", limit=5)
            assert len(results) >= 1
            assert any(result["sysName"] == "rpoA" for result in results)

    def test_find_essential_genes_integration(self):
        """Test find_essential_genes with realistic mock data."""
        mock_genes = {
            "Atu0001": {
                "locusId": "Atu0001",
                "sysName": "essential1",
                "description": "Essential gene for growth",
                "fitness_values": [1.5, 1.2, 0.8],
            },
            "Atu0002": {
                "locusId": "Atu0002",
                "sysName": "nonessential1",
                "description": "Non-essential gene",
                "fitness_values": [0.1, -0.2, 0.0],
            },
        }

        mock_conditions = ["stress_cond", "carbon_cond", "control"]

        def mock_get_gene_fitness(gene_id, condition_filter=None):
            if gene_id not in mock_genes:
                return {"error": "Gene not found"}

            gene_data = mock_genes[gene_id]
            fitness_data = []

            for i, condition in enumerate(mock_conditions):
                if (
                    condition_filter is None
                    or condition_filter.lower() in condition.lower()
                ):
                    if i < len(gene_data["fitness_values"]):
                        fitness_data.append(
                            {
                                "condition": condition,
                                "fitness": gene_data["fitness_values"][i],
                                "description": f"Description for {condition}",
                            }
                        )

            return {
                "gene": {
                    "locusId": gene_data["locusId"],
                    "sysName": gene_data["sysName"],
                    "description": gene_data["description"],
                },
                "fitness_data": fitness_data,
            }

        def mock_interpret_score(score):
            if score is None:
                return {"effect": "unknown", "interpretation": "No data"}
            elif score >= 0.5:
                return {
                    "effect": "gene_benefits_growth",
                    "interpretation": "Essential for growth",
                }
            else:
                return {"effect": "neutral", "interpretation": "Not essential"}

        with (
            patch.object(fitness_loader, "load_data"),
            patch.object(fitness_loader, "genes", mock_genes),
            patch.object(
                fitness_loader, "get_gene_fitness", side_effect=mock_get_gene_fitness
            ),
            patch.object(
                fitness_loader,
                "interpret_fitness_score",
                side_effect=mock_interpret_score,
            ),
        ):
            results = find_essential_genes(
                condition_filter="stress", min_fitness_threshold=0.5, limit=10
            )

            assert isinstance(results, list)
            # Should find Atu0001 as essential
            essential_gene_ids = [result["gene"]["locusId"] for result in results]
            assert "Atu0001" in essential_gene_ids

            # Verify structure
            for result in results:
                assert "gene" in result
                assert "essential_in_conditions" in result
                assert "num_essential_conditions" in result
                assert result["num_essential_conditions"] > 0

    def test_analyze_gene_fitness_integration(self):
        """Test analyze_gene_fitness with comprehensive mock data."""
        mock_fitness_data = {
            "gene": {
                "locusId": "Atu0001",
                "sysName": "testGene",
                "description": "Test gene for analysis",
            },
            "fitness_data": [
                {
                    "condition": "beneficial_cond",
                    "fitness": 0.8,
                },  # Gene benefits growth
                {
                    "condition": "inhibitory_cond",
                    "fitness": -0.8,
                },  # Gene inhibits growth
                {"condition": "neutral_cond1", "fitness": 0.1},  # Neutral
                {"condition": "neutral_cond2", "fitness": -0.1},  # Neutral
                {"condition": "missing_data", "fitness": None},  # No data
            ],
        }

        with patch.object(fitness_loader, "get_gene_fitness") as mock_get:
            mock_get.return_value = mock_fitness_data

            result = analyze_gene_fitness("Atu0001")

            assert "gene" in result
            assert "analysis" in result

            analysis = result["analysis"]
            assert "conditions_where_gene_benefits_growth" in analysis
            assert "conditions_where_gene_inhibits_growth" in analysis
            assert "neutral_conditions" in analysis
            assert "summary" in analysis

            # Verify categorization
            beneficial = analysis["conditions_where_gene_benefits_growth"]
            inhibitory = analysis["conditions_where_gene_inhibits_growth"]
            neutral = analysis["neutral_conditions"]

            assert len(beneficial) == 1
            assert beneficial[0]["fitness"] == 0.8

            assert len(inhibitory) == 1
            assert inhibitory[0]["fitness"] == -0.8

            assert len(neutral) == 2

            # Verify summary
            summary = analysis["summary"]
            assert summary["beneficial_count"] == 1
            assert summary["inhibitory_count"] == 1
            assert summary["neutral_count"] == 2
            assert summary["total_conditions_tested"] == 4  # Excludes None values

    def test_module_integration(self):
        """Test module-related functions integration."""
        mock_gene_to_modules = {
            "Atu0001": [
                {
                    "locus_tag": "Atu0001",
                    "module_id": 1,
                    "gene_weight": 0.8,
                    "product": "Test product",
                    "module_name": "Test Module",
                    "module_category": "Test Category",
                }
            ]
        }

        mock_module_meta = {
            1: {
                "module_id": 1,
                "name": "Test Module",
                "category": "Test Category",
                "count": 5,
            }
        }

        with (
            patch.object(module_loader, "get_modules_for_gene") as mock_get_modules,
            patch.object(module_loader, "search_modules_by_name") as mock_search,
        ):
            mock_get_modules.return_value = mock_gene_to_modules["Atu0001"]
            mock_search.return_value = [
                {
                    "module": mock_module_meta[1],
                    "genes": mock_gene_to_modules["Atu0001"],
                    "gene_count": 1,
                }
            ]

            # Test get_gene_modules
            result = get_gene_modules("Atu0001")
            assert result["gene_id"] == "Atu0001"
            assert len(result["modules"]) == 1
            assert result["module_count"] == 1

            # Test search_modules
            results = search_modules("Test", limit=10)
            assert len(results) == 1
            assert results[0]["module"]["name"] == "Test Module"

    def test_error_handling_integration(self):
        """Test error handling across the integration."""
        # Test with non-existent gene
        result = get_gene_info("NONEXISTENT_GENE")
        assert "error" in result

        # Test gene fitness with non-existent gene
        with patch.object(fitness_loader, "get_gene_fitness") as mock_get:
            mock_get.return_value = {"error": "Gene not found"}

            result = get_gene_fitness("NONEXISTENT_GENE")
            assert "error" in result

        # Test analyze_gene_fitness with error
        with patch.object(fitness_loader, "get_gene_fitness") as mock_get:
            mock_get.return_value = {"error": "Gene not found"}

            result = analyze_gene_fitness("NONEXISTENT_GENE")
            assert "error" in result

    def test_data_consistency_integration(self):
        """Test data consistency across different functions."""
        mock_gene_data = {
            "locusId": "Atu0001",
            "sysName": "testGene",
            "description": "Consistent test gene",
        }

        mock_fitness_response = {
            "gene": mock_gene_data,
            "fitness_data": [{"condition": "test_cond", "fitness": 0.5}],
            "total_conditions": 1,
        }

        with (
            patch.object(fitness_loader, "get_gene_info") as mock_info,
            patch.object(fitness_loader, "get_gene_fitness") as mock_fitness,
        ):
            mock_info.return_value = mock_gene_data
            mock_fitness.return_value = mock_fitness_response

            # Test that gene info is consistent
            info_result = get_gene_info("Atu0001")
            fitness_result = get_gene_fitness("Atu0001")

            assert info_result["locusId"] == fitness_result["gene"]["locusId"]
            assert info_result["sysName"] == fitness_result["gene"]["sysName"]
            assert info_result["description"] == fitness_result["gene"]["description"]
