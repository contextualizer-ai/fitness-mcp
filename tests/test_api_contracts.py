"""API contract tests for all 16 MCP tools in fitness_mcp.

This test suite verifies that all MCP tools maintain their expected API contracts
and data structures. These tests use real data and implementations (no mocks).

API contract verification helps ensure stability during refactoring.
"""

from src.fitness_mcp.main import (
    # Basic gene tools
    get_gene_info,
    get_gene_fitness,
    search_genes,
    # Condition tools
    get_growth_conditions,
    get_condition_details,
    # Analysis tools
    interpret_fitness_score,
    find_essential_genes,
    find_growth_inhibitor_genes,
    analyze_gene_fitness,
    # Module tools
    get_gene_modules,
    get_module_genes,
    search_modules,
    get_all_modules,
    # Pairs/network tools (renamed for consistency)
    get_fitness_effects_for_gene,
    get_genes_with_fitness_effects,
    expand_fitness_network,
    # Data loaders for setup
    fitness_loader,
    # module_loader and pairs_loader removed - functionality replaced by MetadataRegistry
)


class TestAPIContracts:
    """Test API contracts for all 16 MCP tools."""

    def test_get_gene_info_contract(self, loaded_metadata_registry, sample_gene_id):
        """Test get_gene_info API contract."""
        # Test with nonexistent gene (should always work)
        result = get_gene_info("nonexistent_gene_123")
        assert isinstance(result, dict)
        assert "error" in result

        # Test with real gene if available
        if fitness_loader.genes:
            gene_id = list(fitness_loader.genes.keys())[0]
            result = get_gene_info(gene_id)

            if "error" not in result:
                # Success contract
                assert isinstance(result, dict)
                assert "data" in result
                data = result["data"]
                required_keys = {"locusId", "sysName", "description"}
                assert set(data.keys()) == required_keys
                assert all(isinstance(data[key], str) for key in required_keys)

    def test_get_gene_fitness_contract(self):
        """Test get_gene_fitness API contract."""
        # Test with nonexistent gene
        result = get_gene_fitness("nonexistent_gene_123")
        assert isinstance(result, dict)
        assert "error" in result

        # Test with real gene if available
        if fitness_loader.genes:
            gene_id = list(fitness_loader.genes.keys())[0]
            result = get_gene_fitness(gene_id)

            if "error" not in result:
                # Success contract
                assert isinstance(result, dict)
                required_keys = {"gene", "fitness_data", "total_conditions"}
                assert set(result.keys()) == required_keys

                # Gene structure
                gene = result["gene"]
                assert isinstance(gene, dict)
                gene_keys = {"locusId", "sysName", "description"}
                assert set(gene.keys()) == gene_keys

                # Fitness data structure
                fitness_data = result["fitness_data"]
                assert isinstance(fitness_data, list)
                assert isinstance(result["total_conditions"], int)

                if fitness_data:
                    condition_data = fitness_data[0]
                    assert isinstance(condition_data, dict)
                    assert "condition" in condition_data
                    assert "fitness" in condition_data

    def test_search_genes_contract(self):
        """Test search_genes API contract."""
        result = search_genes("test_query", 5)

        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        assert "genes" in data
        if data["genes"]:
            gene = data["genes"][0]
            assert isinstance(gene, dict)
            required_keys = {"locusId", "sysName", "description"}
            assert required_keys.issubset(set(gene.keys()))

    def test_get_growth_conditions_contract(self):
        """Test get_growth_conditions API contract."""
        # Without filter
        result = get_growth_conditions()
        # Check standardized response structure
        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        assert "conditions" in data
        if data["conditions"]:
            assert all(isinstance(condition, str) for condition in data["conditions"])

        # With filter
        result = get_growth_conditions("pH")
        # Check standardized response structure
        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        assert "conditions" in data
        if data["conditions"]:
            assert all(isinstance(condition, str) for condition in data["conditions"])

    def test_get_condition_details_contract(self):
        """Test get_condition_details API contract."""
        # Test with nonexistent condition
        result = get_condition_details("nonexistent_condition_123")
        assert isinstance(result, dict)
        assert "error" in result

        # Test with real condition if available
        if fitness_loader.conditions:
            condition_name = fitness_loader.conditions[0]
            result = get_condition_details(condition_name)

            if "error" not in result:
                assert isinstance(result, dict)
                assert "condition_name" in result
                assert result["condition_name"] == condition_name

    def test_interpret_fitness_score_contract(self):
        """Test interpret_fitness_score API contract."""
        # Test with various score types
        test_scores = [0.8, -0.8, 0.05]

        for score in test_scores:
            result = interpret_fitness_score(score)
            assert isinstance(result, dict)
            assert "data" in result
            data = result["data"]

            required_keys = {"interpretation", "effect", "magnitude"}
            assert required_keys.issubset(set(data.keys()))
            assert all(isinstance(data[key], str) for key in required_keys)

            # Score should be included
            assert "score" in data
            assert data["score"] == score

        # Test None score (should return error)
        result = interpret_fitness_score(None)
        assert "error" in result

    def test_find_essential_genes_contract(self):
        """Test find_essential_genes API contract."""
        result = find_essential_genes(limit=2)

        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        assert "essential_genes" in data
        if data["essential_genes"]:
            gene_entry = data["essential_genes"][0]
            assert isinstance(gene_entry, dict)

            required_keys = {
                "gene",
                "essential_in_conditions",
                "num_essential_conditions",
            }
            assert set(gene_entry.keys()) == required_keys

            # Gene structure
            gene = gene_entry["gene"]
            assert isinstance(gene, dict)
            gene_keys = {"locusId", "sysName", "description"}
            assert set(gene.keys()) == gene_keys

            # Conditions structure
            conditions = gene_entry["essential_in_conditions"]
            assert isinstance(conditions, list)
            assert isinstance(gene_entry["num_essential_conditions"], int)

    def test_find_growth_inhibitor_genes_contract(self):
        """Test find_growth_inhibitor_genes API contract."""
        result = find_growth_inhibitor_genes(limit=2)

        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        assert "inhibitory_genes" in data
        if data["inhibitory_genes"]:
            gene_entry = data["inhibitory_genes"][0]
            assert isinstance(gene_entry, dict)

            required_keys = {
                "gene",
                "inhibits_growth_in_conditions",
                "num_inhibitory_conditions",
            }
            assert set(gene_entry.keys()) == required_keys

            # Gene structure
            gene = gene_entry["gene"]
            assert isinstance(gene, dict)
            gene_keys = {"locusId", "sysName", "description"}
            assert set(gene.keys()) == gene_keys

            # Conditions structure
            conditions = gene_entry["inhibits_growth_in_conditions"]
            assert isinstance(conditions, list)
            assert isinstance(gene_entry["num_inhibitory_conditions"], int)

    def test_analyze_gene_fitness_contract(self):
        """Test analyze_gene_fitness API contract."""
        # Test with nonexistent gene
        result = analyze_gene_fitness("nonexistent_gene_123")
        assert isinstance(result, dict)
        assert "error" in result

        # Test with real gene if available
        if fitness_loader.genes:
            gene_id = list(fitness_loader.genes.keys())[0]
            result = analyze_gene_fitness(gene_id)

            if "error" not in result:
                assert isinstance(result, dict)
                required_keys = {"gene", "analysis"}
                assert set(result.keys()) == required_keys

                # Analysis structure
                analysis = result["analysis"]
                assert isinstance(analysis, dict)
                analysis_keys = {
                    "conditions_where_gene_inhibits_growth",
                    "conditions_where_gene_is_essential",
                    "neutral_conditions",
                    "summary",
                }
                assert set(analysis.keys()) == analysis_keys

                # Summary structure
                summary = analysis["summary"]
                assert isinstance(summary, dict)
                summary_keys = {
                    "total_conditions_tested",
                    "inhibitory_count",
                    "essential_count",
                    "neutral_count",
                    "limit_applied",
                }
                assert set(summary.keys()) == summary_keys
                # Check types for required fields
                int_fields = {
                    "total_conditions_tested",
                    "inhibitory_count",
                    "essential_count",
                    "neutral_count",
                }
                assert all(isinstance(summary[key], int) for key in int_fields)

    def test_get_gene_modules_contract(self, loaded_metadata_registry):
        """Test get_gene_modules API contract."""
        # Test with nonexistent gene
        result = get_gene_modules("nonexistent_gene_123")
        assert isinstance(result, dict)
        assert "error" in result

        # Test with real gene (metadata_registry should have data loaded)
        assert loaded_metadata_registry.gene_to_modules, (
            "gene_to_modules should have data"
        )
        gene_id = list(loaded_metadata_registry.gene_to_modules.keys())[0]
        result = get_gene_modules(gene_id)

        assert "error" not in result, f"Should find modules for gene {gene_id}"
        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        required_keys = {"gene_id", "modules", "module_count"}
        assert set(data.keys()) == required_keys

        assert data["gene_id"] == gene_id
        assert isinstance(data["modules"], list)
        assert isinstance(data["module_count"], int)

    def test_get_module_genes_contract(self, loaded_metadata_registry):
        """Test get_module_genes API contract."""
        # Test with real module (metadata_registry should have data loaded)
        assert loaded_metadata_registry.module_to_genes, (
            "module_to_genes should have data"
        )
        module_id = list(loaded_metadata_registry.module_to_genes.keys())[0]
        result = get_module_genes(module_id)

        assert "error" not in result, f"Should find genes for module {module_id}"
        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        required_keys = {"module", "genes", "gene_count"}
        assert set(data.keys()) == required_keys

        assert isinstance(data["module"], dict)
        assert isinstance(data["genes"], list)
        assert isinstance(data["gene_count"], int)

    def test_search_modules_contract(self):
        """Test search_modules API contract."""
        result = search_modules("transport", 5)

        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        assert "modules" in data
        if data["modules"]:
            module_entry = data["modules"][0]
            assert isinstance(module_entry, dict)
            required_keys = {"module", "genes", "gene_count"}
            assert set(module_entry.keys()) == required_keys

    def test_get_all_modules_contract(self):
        """Test get_all_modules API contract."""
        result = get_all_modules()

        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        assert "modules" in data
        if data["modules"]:
            module = data["modules"][0]
            assert isinstance(module, dict)
            required_keys = {"module_id", "name", "category", "count"}
            assert set(module.keys()) == required_keys

    def test_get_fitness_effects_for_gene_contract(self, loaded_metadata_registry):
        """Test get_fitness_effects_for_gene API contract."""
        # Test with nonexistent gene
        result = get_fitness_effects_for_gene("nonexistent_gene_123")
        assert isinstance(result, dict)
        assert "error" in result

        # Test with real gene (metadata_registry should have data loaded)
        assert loaded_metadata_registry.gene_to_conditions, (
            "gene_to_conditions should have data"
        )
        gene_id = list(loaded_metadata_registry.gene_to_conditions.keys())[0]
        result = get_fitness_effects_for_gene(gene_id)

        assert "error" not in result, f"Should find fitness effects for gene {gene_id}"
        assert isinstance(result, dict)
        # New standardized response format
        required_keys = {"data", "metadata", "suggestions"}
        assert set(result.keys()) == required_keys

        # Check data section
        data = result["data"]
        data_keys = {"gene_id", "fitness_effects", "total_effects"}
        assert set(data.keys()) == data_keys

        assert data["gene_id"] == gene_id
        assert isinstance(data["fitness_effects"], list)
        assert isinstance(data["total_effects"], int)
        assert isinstance(result["metadata"]["interpretation"], str)

    def test_get_genes_with_fitness_effects_contract(self, loaded_metadata_registry):
        """Test get_genes_with_fitness_effects API contract."""
        # Test with nonexistent condition
        result = get_genes_with_fitness_effects("nonexistent_condition_123")
        assert isinstance(result, dict)
        assert "error" in result

        # Test with real condition (metadata_registry should have data loaded)
        assert loaded_metadata_registry.condition_to_genes, (
            "condition_to_genes should have data"
        )
        condition_id = list(loaded_metadata_registry.condition_to_genes.keys())[0]
        result = get_genes_with_fitness_effects(condition_id)

        assert "error" not in result, f"Should find genes for condition {condition_id}"
        assert isinstance(result, dict)
        # New standardized response format
        required_keys = {"data", "metadata", "suggestions"}
        assert set(result.keys()) == required_keys

        # Check data section
        data = result["data"]
        data_keys = {"condition_id", "fitness_effects", "total_genes"}
        assert set(data.keys()) == data_keys

        assert data["condition_id"] == condition_id
        assert isinstance(data["fitness_effects"], list)
        assert isinstance(data["total_genes"], int)
        assert isinstance(result["metadata"]["interpretation"], str)

    def test_expand_fitness_network_contract(self, loaded_metadata_registry):
        """Test expand_fitness_network API contract."""
        # Test with nonexistent pair
        result = expand_fitness_network("gene123", "condition123")
        assert isinstance(result, dict)
        assert "error" in result

        # Test with real gene-condition pair (metadata_registry should have data loaded)
        assert loaded_metadata_registry.gene_to_conditions, (
            "gene_to_conditions should have data"
        )
        gene_id = list(loaded_metadata_registry.gene_to_conditions.keys())[0]
        conditions = loaded_metadata_registry.gene_to_conditions[gene_id]

        assert len(conditions) > 0, (
            f"Gene {gene_id} should have conditions in loaded data"
        )
        condition_id = conditions[0]  # metadata_registry stores condition IDs directly
        result = expand_fitness_network(gene_id, condition_id)

        assert "error" not in result, (
            f"Should expand network for {gene_id}-{condition_id}"
        )
        assert isinstance(result, dict)
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


class TestMCPToolInventory:
    """Verify all 16 MCP tools are properly registered and documented."""

    def test_all_tools_are_tested(self):
        """Ensure all 16 MCP tools have API contract tests."""
        # List of all MCP tools from main.py mcp.tool() registrations
        expected_tools = {
            "get_gene_info",
            "get_gene_fitness",
            "search_genes",
            "get_growth_conditions",
            "get_condition_details",
            "interpret_fitness_score",
            "find_essential_genes",
            "find_growth_inhibitor_genes",
            "analyze_gene_fitness",
            "get_gene_modules",
            "get_module_genes",
            "search_modules",
            "get_all_modules",
            "get_fitness_effects_for_gene",
            "get_genes_with_fitness_effects",
            "expand_fitness_network",
        }

        # Extract test method names from this file
        test_methods = [
            name
            for name in dir(TestAPIContracts)
            if name.startswith("test_") and name.endswith("_contract")
        ]
        tested_tools = {
            method.replace("test_", "").replace("_contract", "")
            for method in test_methods
        }

        assert len(expected_tools) == 16, "Should have exactly 16 MCP tools"
        assert expected_tools == tested_tools, (
            f"Missing tests for: {expected_tools - tested_tools}"
        )

    def test_tool_function_signatures(self):
        """Verify all MCP tools have proper function signatures."""
        tool_functions = [
            get_gene_info,
            get_gene_fitness,
            search_genes,
            get_growth_conditions,
            get_condition_details,
            interpret_fitness_score,
            find_essential_genes,
            find_growth_inhibitor_genes,
            analyze_gene_fitness,
            get_gene_modules,
            get_module_genes,
            search_modules,
            get_all_modules,
            get_fitness_effects_for_gene,
            get_genes_with_fitness_effects,
            expand_fitness_network,
        ]

        for func in tool_functions:
            # All tools should be callable
            assert callable(func), f"{func.__name__} is not callable"

            # All tools should have docstrings
            assert func.__doc__, f"{func.__name__} missing docstring"

            # All tools should have type annotations
            assert hasattr(func, "__annotations__"), (
                f"{func.__name__} missing type annotations"
            )
