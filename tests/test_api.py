"""Basic API tests for fitness_mcp.

Tests use real data and session fixtures for performance.
No mocks or patches are used per project policy.
"""

from src.fitness_mcp.main import (
    get_gene_info,
    get_gene_fitness,
    search_genes,
    get_growth_conditions,
    interpret_fitness_score,
)


def test_reality():
    """Basic sanity check."""
    assert 1 == 1


def test_gene_info_api_basic(loaded_metadata_registry):
    """Test basic gene info API functionality with real data."""
    # Test with a real gene from the loaded data
    if loaded_metadata_registry.genes:
        gene_id = list(loaded_metadata_registry.genes.keys())[0]
        result = get_gene_info(gene_id)

        if "error" not in result:
            assert "locusId" in result
            assert "sysName" in result
            assert "description" in result
            assert isinstance(result["locusId"], str)
            assert isinstance(result["sysName"], str)
            assert isinstance(result["description"], str)

    # Test with non-existent gene
    result = get_gene_info("nonexistent_gene_123")
    assert "error" in result


def test_gene_fitness_api_basic(loaded_metadata_registry):
    """Test basic gene fitness API functionality with real data."""
    # Test with a real gene from the loaded data
    if loaded_metadata_registry.genes:
        gene_id = list(loaded_metadata_registry.genes.keys())[0]
        result = get_gene_fitness(gene_id)

        if "error" not in result:
            assert "gene" in result
            assert "fitness_data" in result
            assert "total_conditions" in result
            assert isinstance(result["fitness_data"], list)

            # Check gene structure
            gene = result["gene"]
            assert "locusId" in gene
            assert "sysName" in gene
            assert "description" in gene

    # Test with non-existent gene
    result = get_gene_fitness("nonexistent_gene_123")
    assert "error" in result


def test_search_genes_api_basic(loaded_metadata_registry):
    """Test basic search genes API functionality with real data."""
    # Data is already loaded via fixture
    assert loaded_metadata_registry.loaded, "Fitness data must be loaded for API tests"
    assert len(loaded_metadata_registry.genes) > 0, (
        "Fitness data must contain genes for API tests"
    )

    result = search_genes("gene")  # Search for common term

    assert isinstance(result, list)
    # Test structure if results exist
    if len(result) > 0:
        assert all("locusId" in gene for gene in result)
        assert all("sysName" in gene for gene in result)
        assert all("description" in gene for gene in result)


def test_growth_conditions_api_basic(loaded_metadata_registry):
    """Test basic growth conditions API functionality with real data."""
    # Data is already loaded via fixture
    assert loaded_metadata_registry.loaded, "Fitness data must be loaded for API tests"
    assert len(loaded_metadata_registry.conditions) > 0, (
        "Fitness data must contain conditions for API tests"
    )

    result = get_growth_conditions()

    assert isinstance(result, list)
    assert len(result) > 0  # Should have conditions with real data
    assert all(isinstance(condition, str) for condition in result)


def test_interpret_fitness_score_api_basic():
    """Test basic fitness score interpretation API with real data."""
    # Test positive score (gene inhibits growth)
    result = interpret_fitness_score(0.8)
    assert "interpretation" in result
    assert "effect" in result
    assert "magnitude" in result
    assert "score" in result
    assert result["score"] == 0.8
    assert result["effect"] == "gene_inhibits_growth"

    # Test negative score (gene is essential)
    result = interpret_fitness_score(-0.8)
    assert result["effect"] == "gene_benefits_growth"

    # Test neutral score
    result = interpret_fitness_score(0.05)
    assert result["effect"] == "neutral"

    # Test None score
    result = interpret_fitness_score(None)
    assert result["effect"] == "unknown"
