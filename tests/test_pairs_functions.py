"""Tests for pairs/fitness effects functions in fitness_mcp.main module.

Tests use real data and session fixtures for performance.
No mocks or patches are used per project policy.
"""

from src.fitness_mcp.main import (
    get_fitness_effects_for_gene,
    get_genes_with_fitness_effects,
    expand_fitness_network,
)


class TestPairsFunctions:
    """Test cases for pairs/fitness effects MCP tool functions."""

    def test_get_fitness_effects_for_gene_with_real_data(
        self, loaded_metadata_registry
    ):
        """Test get_fitness_effects_for_gene with real data."""
        # Data is already loaded via fixture - must have gene-condition data
        assert len(loaded_metadata_registry.gene_to_conditions) > 0, (
            "Must have gene-condition fitness effects loaded"
        )

        # Use a real gene from the loaded data instead of hardcoded one
        gene_id = list(loaded_metadata_registry.gene_to_conditions.keys())[0]
        result = get_fitness_effects_for_gene(gene_id)

        # With real data, this should succeed
        assert "error" not in result, (
            f"Expected successful result for gene {gene_id}, got: {result}"
        )

        # Check standardized response structure
        assert "data" in result
        assert "metadata" in result
        assert "suggestions" in result

        data = result["data"]
        assert "gene_id" in data
        assert "fitness_effects" in data
        assert "total_effects" in data

        assert data["gene_id"] == gene_id
        assert isinstance(data["fitness_effects"], list)
        assert isinstance(data["total_effects"], int)
        assert len(data["fitness_effects"]) > 0, (
            "Should have fitness effects for genes with significant effects"
        )

        # Check structure of first fitness effect
        effect = data["fitness_effects"][0]
        assert "condition" in effect
        assert "fitness_value" in effect
        assert isinstance(effect["fitness_value"], (int, float))

    def test_get_fitness_effects_for_gene_nonexistent(self, loaded_metadata_registry):
        """Test get_fitness_effects_for_gene with non-existent gene."""
        # Data is already loaded via fixture
        result = get_fitness_effects_for_gene("NonExistentGene123")

        assert "error" in result
        assert "NonExistentGene123" in result["error"]

    def test_get_genes_with_fitness_effects_with_real_data(
        self, loaded_metadata_registry
    ):
        """Test get_genes_with_fitness_effects with real data."""
        # Data is already loaded via fixture - must have condition-gene data
        assert len(loaded_metadata_registry.condition_to_genes) > 0, (
            "Must have condition-gene fitness effects loaded"
        )

        # Get a real condition from the data
        condition_id = list(loaded_metadata_registry.condition_to_genes.keys())[0]
        result = get_genes_with_fitness_effects(condition_id)

        # With real data, this should succeed
        assert "error" not in result, (
            f"Expected successful result for condition {condition_id}, got: {result}"
        )

        # Check standardized response structure
        assert "data" in result
        assert "metadata" in result
        assert "suggestions" in result

        data = result["data"]
        assert "condition_id" in data
        assert "fitness_effects" in data
        assert "total_genes" in data

        assert data["condition_id"] == condition_id
        assert isinstance(data["fitness_effects"], list)
        assert isinstance(data["total_genes"], int)
        assert len(data["fitness_effects"]) > 0, (
            "Should have genes for conditions with fitness effects"
        )

        # Check structure of first gene effect
        gene_effect = data["fitness_effects"][0]
        assert "gene" in gene_effect
        assert "fitness_value" in gene_effect
        assert isinstance(gene_effect["fitness_value"], (int, float))

    def test_get_genes_with_fitness_effects_nonexistent(self, loaded_metadata_registry):
        """Test get_genes_with_fitness_effects with non-existent condition."""
        # Data is already loaded via fixture
        result = get_genes_with_fitness_effects("NonExistentCondition123")

        assert "error" in result
        assert "NonExistentCondition123" in result["error"]

    def test_expand_fitness_network_with_real_data(self, loaded_metadata_registry):
        """Test expand_fitness_network with real data."""
        # Data is already loaded via fixture - must have gene-condition pairs
        assert len(loaded_metadata_registry.gene_to_conditions) > 0, (
            "Must have gene-condition fitness effects loaded"
        )

        # Find a real gene-condition pair
        gene_id = list(loaded_metadata_registry.gene_to_conditions.keys())[0]
        conditions = loaded_metadata_registry.gene_to_conditions[gene_id]
        assert len(conditions) > 0, "Gene must have associated conditions"

        condition_id = conditions[0]
        result = expand_fitness_network(gene_id, condition_id)

        # With real data, this should succeed
        assert "error" not in result, (
            f"Expected successful result for {gene_id}-{condition_id}, got: {result}"
        )

        # Verify structure
        assert "query" in result
        assert "first_hop" in result
        assert "second_hop" in result
        assert "network_size" in result
        assert "interpretation" in result

        # Verify query details
        query = result["query"]
        assert query["gene_id"] == gene_id
        assert query["condition_id"] == condition_id
        assert "fitness_value" in query

        # Verify first hop
        first_hop = result["first_hop"]
        assert "conditions_for_query_gene" in first_hop
        assert "genes_for_query_condition" in first_hop
        assert "num_conditions" in first_hop
        assert "num_genes" in first_hop
        assert first_hop["num_conditions"] > 0, "Should have conditions in network"
        assert first_hop["num_genes"] > 0, "Should have genes in network"

        # Verify second hop
        second_hop = result["second_hop"]
        assert "all_genes_in_network" in second_hop
        assert "all_conditions_in_network" in second_hop
        assert "num_total_genes" in second_hop
        assert "num_total_conditions" in second_hop

        # Verify network size
        network_size = result["network_size"]
        assert "gene_expansion_factor" in network_size
        assert "condition_expansion_factor" in network_size

    def test_expand_fitness_network_nonexistent_pair(self, loaded_metadata_registry):
        """Test expand_fitness_network with non-existent gene-condition pair."""
        # Data is already loaded via fixture
        result = expand_fitness_network("Atu3150", "NonExistentCondition123")

        assert "error" in result
        assert "No significant fitness effect found" in result["error"]


class TestPairsFunctionsBehaviorPreservation:
    """Test cases to preserve exact behavior before refactoring."""

    def test_get_fitness_effects_for_gene_output_schema(self, loaded_metadata_registry):
        """Test that get_fitness_effects_for_gene maintains exact output schema."""
        # Data is already loaded via fixture - must have data
        assert len(loaded_metadata_registry.gene_to_conditions) > 0, (
            "Must have gene-condition fitness effects loaded"
        )

        # Use real data to test schema - should have successful results with real data
        gene_id = list(loaded_metadata_registry.gene_to_conditions.keys())[0]
        result = get_fitness_effects_for_gene(gene_id)

        # With real data and a valid gene_id, we should get successful results
        assert "error" not in result, (
            f"Expected successful result for gene {gene_id}, got: {result}"
        )

        # Standardized response schema
        required_keys = {"data", "metadata", "suggestions"}
        assert set(result.keys()) == required_keys

        # Data section
        data = result["data"]
        data_keys = {"gene_id", "fitness_effects", "total_effects"}
        assert set(data.keys()) == data_keys
        assert isinstance(data["gene_id"], str)
        assert isinstance(data["fitness_effects"], list)
        assert isinstance(data["total_effects"], int)
        assert len(data["fitness_effects"]) > 0, (
            "Should have fitness effects for genes with significant effects"
        )

        # Fitness effects structure
        effect = data["fitness_effects"][0]
        assert "condition" in effect
        assert "fitness_value" in effect

    def test_get_genes_with_fitness_effects_output_schema(
        self, loaded_metadata_registry
    ):
        """Test that get_genes_with_fitness_effects maintains exact output schema."""
        # Data is already loaded via fixture - must have data
        assert len(loaded_metadata_registry.condition_to_genes) > 0, (
            "Must have condition-gene fitness effects loaded"
        )

        # Use real data to test schema - should have successful results with real data
        condition_id = list(loaded_metadata_registry.condition_to_genes.keys())[0]
        result = get_genes_with_fitness_effects(condition_id)

        # With real data and a valid condition_id, we should get successful results
        assert "error" not in result, (
            f"Expected successful result for condition {condition_id}, got: {result}"
        )

        # Standardized response schema
        required_keys = {"data", "metadata", "suggestions"}
        assert set(result.keys()) == required_keys

        # Data section
        data = result["data"]
        data_keys = {"condition_id", "fitness_effects", "total_genes"}
        assert set(data.keys()) == data_keys
        assert isinstance(data["condition_id"], str)
        assert isinstance(data["fitness_effects"], list)
        assert isinstance(data["total_genes"], int)
        assert len(data["fitness_effects"]) > 0, (
            "Should have genes for conditions with fitness effects"
        )

    def test_expand_fitness_network_output_schema(self, loaded_metadata_registry):
        """Test that expand_fitness_network maintains exact output schema."""
        # Data is already loaded via fixture - must have data
        assert len(loaded_metadata_registry.gene_to_conditions) > 0, (
            "Must have gene-condition fitness effects loaded"
        )

        # Use real data to test schema - should have successful results with real data
        gene_id = list(loaded_metadata_registry.gene_to_conditions.keys())[0]
        conditions = loaded_metadata_registry.gene_to_conditions[gene_id]
        assert len(conditions) > 0, "Gene must have associated conditions"

        condition_id = conditions[0]
        result = expand_fitness_network(gene_id, condition_id)

        # With real data and valid gene-condition pair, we should get successful results
        assert "error" not in result, (
            f"Expected successful result for {gene_id}-{condition_id}, got: {result}"
        )

        # Exact schema preservation
        required_keys = {
            "query",
            "first_hop",
            "second_hop",
            "network_size",
            "interpretation",
        }
        assert set(result.keys()) == required_keys

        # Query structure
        query = result["query"]
        query_keys = {"gene_id", "condition_id", "fitness_value"}
        assert set(query.keys()) == query_keys

        # First hop structure
        first_hop = result["first_hop"]
        first_hop_keys = {
            "conditions_for_query_gene",
            "genes_for_query_condition",
            "num_conditions",
            "num_genes",
        }
        assert set(first_hop.keys()) == first_hop_keys

        # Second hop structure
        second_hop = result["second_hop"]
        second_hop_keys = {
            "all_genes_in_network",
            "all_conditions_in_network",
            "num_total_genes",
            "num_total_conditions",
        }
        assert set(second_hop.keys()) == second_hop_keys

        # Network size structure
        network_size = result["network_size"]
        network_size_keys = {
            "gene_expansion_factor",
            "condition_expansion_factor",
        }
        assert set(network_size.keys()) == network_size_keys
