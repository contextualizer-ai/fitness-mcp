"""Tests for pairs/fitness effects functions in fitness_mcp.main module.

Tests use real data loaders and real implementations.
No mocks or patches are used per project policy.
"""

from src.fitness_mcp.main import (
    get_conditions_for_gene,
    get_genes_for_condition,
    expand_gene_condition_network,
    pairs_loader,
)


class TestPairsFunctions:
    """Test cases for pairs/fitness effects MCP tool functions."""

    def test_get_conditions_for_gene_with_real_data(self):
        """Test get_conditions_for_gene with real data loader."""
        pairs_loader.load_data()

        # Test with a gene that should have conditions
        result = get_conditions_for_gene("Atu3150")

        # Check basic structure
        assert "gene_id" in result
        assert "conditions" in result
        assert "total_conditions" in result
        assert "interpretation" in result

        # If the gene exists, verify data types
        if "error" not in result:
            assert result["gene_id"] == "Atu3150"
            assert isinstance(result["conditions"], list)
            assert isinstance(result["total_conditions"], int)
            assert isinstance(result["interpretation"], str)

            # If conditions exist, check structure
            if result["conditions"]:
                condition = result["conditions"][0]
                assert "condition" in condition
                assert "value" in condition
                assert isinstance(condition["value"], (int, float))

    def test_get_conditions_for_gene_nonexistent(self):
        """Test get_conditions_for_gene with non-existent gene."""
        pairs_loader.load_data()

        result = get_conditions_for_gene("NonExistentGene123")

        assert "error" in result
        assert "NonExistentGene123" in result["error"]

    def test_get_genes_for_condition_with_real_data(self):
        """Test get_genes_for_condition with real data loader."""
        pairs_loader.load_data()

        # Get a real condition from the data
        if pairs_loader.condition_to_genes:
            condition_id = list(pairs_loader.condition_to_genes.keys())[0]

            result = get_genes_for_condition(condition_id)

            # Check basic structure
            assert "condition_id" in result
            assert "genes" in result
            assert "total_genes" in result
            assert "interpretation" in result

            if "error" not in result:
                assert result["condition_id"] == condition_id
                assert isinstance(result["genes"], list)
                assert isinstance(result["total_genes"], int)
                assert isinstance(result["interpretation"], str)

                # If genes exist, check structure
                if result["genes"]:
                    gene = result["genes"][0]
                    assert "gene" in gene
                    assert "value" in gene
                    assert isinstance(gene["value"], (int, float))

    def test_get_genes_for_condition_nonexistent(self):
        """Test get_genes_for_condition with non-existent condition."""
        pairs_loader.load_data()

        result = get_genes_for_condition("NonExistentCondition123")

        assert "error" in result
        assert "NonExistentCondition123" in result["error"]

    def test_expand_gene_condition_network_with_real_data(self):
        """Test expand_gene_condition_network with real data."""
        pairs_loader.load_data()

        # Find a real gene-condition pair
        if pairs_loader.gene_to_conditions:
            gene_id = list(pairs_loader.gene_to_conditions.keys())[0]
            conditions = pairs_loader.gene_to_conditions[gene_id]

            if conditions:
                condition_id = conditions[0]["condition"]

                result = expand_gene_condition_network(gene_id, condition_id)

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

    def test_expand_gene_condition_network_nonexistent_pair(self):
        """Test expand_gene_condition_network with non-existent gene-condition pair."""
        pairs_loader.load_data()

        result = expand_gene_condition_network("Atu3150", "NonExistentCondition123")

        assert "error" in result
        assert "No significant fitness value found" in result["error"]


class TestPairsFunctionsBehaviorPreservation:
    """Test cases to preserve exact behavior before refactoring."""

    def test_get_conditions_for_gene_output_schema(self):
        """Test that get_conditions_for_gene maintains exact output schema."""
        pairs_loader.load_data()

        # Use real data to test schema
        if pairs_loader.gene_to_conditions:
            gene_id = list(pairs_loader.gene_to_conditions.keys())[0]
            result = get_conditions_for_gene(gene_id)

            # Exact schema preservation
            if "error" not in result:
                required_keys = {
                    "gene_id",
                    "conditions",
                    "total_conditions",
                    "interpretation",
                }
                assert set(result.keys()) == required_keys

                # Value types
                assert isinstance(result["gene_id"], str)
                assert isinstance(result["conditions"], list)
                assert isinstance(result["total_conditions"], int)
                assert isinstance(result["interpretation"], str)

                # Conditions structure
                if result["conditions"]:
                    condition = result["conditions"][0]
                    assert "condition" in condition
                    assert "value" in condition

    def test_get_genes_for_condition_output_schema(self):
        """Test that get_genes_for_condition maintains exact output schema."""
        pairs_loader.load_data()

        # Use real data to test schema
        if pairs_loader.condition_to_genes:
            condition_id = list(pairs_loader.condition_to_genes.keys())[0]
            result = get_genes_for_condition(condition_id)

            # Exact schema preservation
            if "error" not in result:
                required_keys = {
                    "condition_id",
                    "genes",
                    "total_genes",
                    "interpretation",
                }
                assert set(result.keys()) == required_keys

                # Value types
                assert isinstance(result["condition_id"], str)
                assert isinstance(result["genes"], list)
                assert isinstance(result["total_genes"], int)
                assert isinstance(result["interpretation"], str)

    def test_expand_gene_condition_network_output_schema(self):
        """Test that expand_gene_condition_network maintains exact output schema."""
        pairs_loader.load_data()

        # Use real data to test schema
        if pairs_loader.gene_to_conditions:
            gene_id = list(pairs_loader.gene_to_conditions.keys())[0]
            conditions = pairs_loader.gene_to_conditions[gene_id]

            if conditions:
                condition_id = conditions[0]["condition"]
                result = expand_gene_condition_network(gene_id, condition_id)

                # Exact schema preservation
                if "error" not in result:
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
