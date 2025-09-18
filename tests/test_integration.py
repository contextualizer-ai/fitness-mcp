"""Integration tests for fitness_mcp using real data files."""

import pytest
from unittest.mock import patch

from src.fitness_mcp.main import (
    fitness_loader,
    # module_loader removed - functionality replaced by MetadataRegistry
    metadata_registry,
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
    def setup_integration_environment(self, loaded_metadata_registry):
        """Set up environment for integration tests using shared fixture."""
        # The loaded_metadata_registry fixture ensures data is loaded
        # and will skip tests if data is not available
        yield

    def test_fitness_loader_real_data(self, loaded_metadata_registry):
        """Test fitness loader with real data files."""
        # Data is already loaded via fixture
        fitness_file = "data/fit_t.tab"
        exp_file = "data/exp_organism_Agro.txt"
        # Temporarily override file paths
        original_data_file = fitness_loader.data_file
        original_exp_file = fitness_loader.exp_desc_file

        try:
            fitness_loader.data_file = fitness_file
            fitness_loader.exp_desc_file = exp_file
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

    # test_module_loader_real_data removed - ModuleDataLoader was replaced by MetadataRegistry
    # Equivalent functionality is tested in tests/test_metadata_registry.py

    def test_search_genes_integration(self):
        """Test search_genes integration with real data."""
        # Test search functionality with real data
        first_gene = next(iter(metadata_registry.genes.values()))

        # Search by gene ID
        results = search_genes(first_gene.locus_id, limit=5)
        assert len(results) >= 1
        assert any(result["locusId"] == first_gene.locus_id for result in results)

        # Search by partial description if available
        if first_gene.description and len(first_gene.description) > 3:
            search_term = first_gene.description[:4].lower()
            results = search_genes(search_term, limit=5)
            assert isinstance(results, list)  # Should return valid results

    def test_find_essential_genes_integration(self):
        """Test find_essential_genes integration."""
        # Test basic functionality - should return a list regardless of data availability
        results = find_essential_genes(
            condition_filter=None, min_fitness_threshold=0.5, limit=10
        )
        assert isinstance(results, list)

    def test_analyze_gene_fitness_integration(self):
        """Test analyze_gene_fitness integration with real data."""
        # Test with a gene that should exist in real data
        first_gene = next(iter(metadata_registry.genes.values()))

        result = analyze_gene_fitness(first_gene.locus_id)

        # Should return proper structure regardless of data content
        assert "gene" in result or "error" in result

        if "gene" in result:
            assert "analysis" in result
            analysis = result["analysis"]
            assert "conditions_where_gene_is_essential" in analysis
            assert "conditions_where_gene_inhibits_growth" in analysis
            assert "neutral_conditions" in analysis
            assert "summary" in analysis

    def test_module_integration(self):
        """Test module-related functions integration."""
        # Test basic functionality - should handle missing data gracefully
        result = get_gene_modules("NONEXISTENT_GENE")
        assert "error" in result or "gene_id" in result

        # Test search_modules
        results = search_modules("test", limit=10)
        assert isinstance(results, list)

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
        # Test that functions handle missing data consistently
        info_result = get_gene_info("NONEXISTENT_GENE")
        fitness_result = get_gene_fitness("NONEXISTENT_GENE")

        # Both should return error when gene doesn't exist
        assert "error" in info_result
        assert "error" in fitness_result
