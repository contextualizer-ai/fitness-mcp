"""Tests for MCP tool functions in fitness_mcp.main module.

Tests use real data and session fixtures for performance.
No mocks or patches are used per project policy.
"""

from src.fitness_mcp.main import (
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
)


class TestMCPToolFunctions:
    """Test cases for MCP tool functions using real data."""

    def test_get_gene_info_success(self, loaded_metadata_registry):
        """Test get_gene_info with existing gene."""
        # Data is already loaded via fixture
        # Test with a known gene if it exists
        if loaded_metadata_registry.genes:
            gene_id = list(loaded_metadata_registry.genes.keys())[0]
            result = get_gene_info(gene_id)

            if "error" not in result:
                # Check standardized response structure
                assert "data" in result
                data = result["data"]
                assert "locusId" in data
                assert "sysName" in data
                assert "description" in data
                assert isinstance(data["locusId"], str)
                assert isinstance(data["sysName"], str)
                assert isinstance(data["description"], str)

    def test_get_gene_info_not_found(self):
        """Test get_gene_info with non-existent gene."""
        result = get_gene_info("nonexistent_gene_123")

        assert "error" in result
        assert "not found" in result["error"]

    def test_get_gene_fitness(self, loaded_metadata_registry):
        """Test get_gene_fitness function."""
        # Data is already loaded via fixture
        # Test with a known gene if it exists
        if loaded_metadata_registry.genes:
            gene_id = list(loaded_metadata_registry.genes.keys())[0]
            result = get_gene_fitness(gene_id)

            if "error" not in result:
                assert "gene" in result
                assert "fitness_data" in result
                assert "total_conditions" in result

                # Check gene structure
                gene = result["gene"]
                assert "locusId" in gene
                assert "sysName" in gene
                assert "description" in gene

                # Check fitness data structure
                fitness_data = result["fitness_data"]
                assert isinstance(fitness_data, list)
                if fitness_data:
                    condition_data = fitness_data[0]
                    assert "condition" in condition_data
                    assert "fitness" in condition_data

    def test_search_genes(self, loaded_metadata_registry):
        """Test search_genes function."""
        # Data is already loaded via fixture
        # Use a common term likely to find results
        result = search_genes("ribosom", 5)

        # Check standardized response structure
        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        assert "genes" in data
        if data["genes"]:
            gene = data["genes"][0]
            assert "locusId" in gene
            assert "sysName" in gene
            assert "description" in gene

    def test_get_growth_conditions(self, loaded_metadata_registry):
        """Test get_growth_conditions function."""
        # Data is already loaded via fixture
        result = get_growth_conditions()

        # Check standardized response structure
        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        assert "conditions" in data
        if data["conditions"]:
            assert isinstance(data["conditions"][0], str)

    def test_get_growth_conditions_with_filter(self, loaded_metadata_registry):
        """Test get_growth_conditions with filter."""
        # Data is already loaded via fixture
        result = get_growth_conditions("pH")

        # Check standardized response structure
        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        assert "conditions" in data
        # All returned conditions should contain the filter term
        for condition in data["conditions"]:
            assert "pH" in condition or "ph" in condition.lower()

    def test_get_condition_details(self, loaded_metadata_registry):
        """Test get_condition_details function."""
        # Data is already loaded via fixture
        # Use a real condition if available
        if loaded_metadata_registry.conditions:
            condition_name = list(loaded_metadata_registry.conditions.keys())[0]
            result = get_condition_details(condition_name)

            if "error" not in result:
                # Check standardized response structure
                assert "data" in result
                data = result["data"]
                assert "condition_info" in data
                assert data["condition_info"]["condition_name"] == condition_name

    def test_interpret_fitness_score(self):
        """Test interpret_fitness_score function."""
        # Test positive score (gene inhibits growth - knockout improves fitness)
        result = interpret_fitness_score(0.8)

        # Check standardized response structure
        assert "data" in result
        data = result["data"]
        assert "interpretation" in data
        assert "effect" in data
        assert "magnitude" in data
        assert "score" in data
        assert data["score"] == 0.8
        assert data["effect"] == "gene_inhibits_growth"

        # Test negative score (gene benefits growth - knockout reduces fitness)
        result = interpret_fitness_score(-0.8)
        assert result["data"]["effect"] == "gene_benefits_growth"

        # Test neutral score
        result = interpret_fitness_score(0.05)
        assert result["data"]["effect"] == "neutral"

        # Test None score
        result = interpret_fitness_score(None)
        # This should return an error for None input
        assert "error" in result

    def test_find_essential_genes(self, loaded_metadata_registry):
        """Test find_essential_genes function."""
        # Data is already loaded via fixture
        result = find_essential_genes(limit=2)

        # Check standardized response structure
        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        assert "essential_genes" in data
        # Result could be empty if no essential genes found with default threshold
        if data["essential_genes"]:
            gene_entry = data["essential_genes"][0]
            assert "gene" in gene_entry
            assert "essential_in_conditions" in gene_entry
            assert "num_essential_conditions" in gene_entry

            # Check gene structure
            gene = gene_entry["gene"]
            assert "locusId" in gene
            assert "sysName" in gene
            assert "description" in gene

    def test_find_growth_inhibitor_genes(self, loaded_metadata_registry):
        """Test find_growth_inhibitor_genes function."""
        # Data is already loaded via fixture
        result = find_growth_inhibitor_genes(limit=2)

        # Check standardized response structure
        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        assert "inhibitory_genes" in data
        # Result could be empty if no inhibitor genes found with default threshold
        if data["inhibitory_genes"]:
            gene_entry = data["inhibitory_genes"][0]
            assert "gene" in gene_entry
            assert "inhibits_growth_in_conditions" in gene_entry
            assert "num_inhibitory_conditions" in gene_entry

    def test_analyze_gene_fitness(self, loaded_metadata_registry):
        """Test analyze_gene_fitness function."""
        # Data is already loaded via fixture
        # Test with a known gene if available
        if loaded_metadata_registry.genes:
            gene_id = list(loaded_metadata_registry.genes.keys())[0]
            result = analyze_gene_fitness(gene_id)

            if "error" not in result:
                assert "gene" in result
                assert "analysis" in result

                analysis = result["analysis"]
                assert "conditions_where_gene_inhibits_growth" in analysis
                assert "conditions_where_gene_is_essential" in analysis
                assert "neutral_conditions" in analysis
                assert "summary" in analysis

                # Check summary structure
                summary = analysis["summary"]
                assert "total_conditions_tested" in summary
                assert "inhibitory_count" in summary
                assert "essential_count" in summary
                assert "neutral_count" in summary

    def test_analyze_gene_fitness_error(self):
        """Test analyze_gene_fitness with error."""
        result = analyze_gene_fitness("nonexistent_gene_123")

        assert "error" in result

    def test_get_gene_modules_success(self, loaded_metadata_registry):
        """Test get_gene_modules with existing gene."""
        # Data is already loaded via fixture
        # Test with a known gene if modules exist
        if loaded_metadata_registry.gene_to_modules:
            gene_id = list(loaded_metadata_registry.gene_to_modules.keys())[0]
            result = get_gene_modules(gene_id)

            if "error" not in result:
                # Check standardized response structure
                assert "data" in result
                data = result["data"]
                assert "gene_id" in data
                assert "modules" in data
                assert "module_count" in data
                assert data["gene_id"] == gene_id
                assert isinstance(data["modules"], list)
                assert isinstance(data["module_count"], int)

    def test_get_gene_modules_not_found(self):
        """Test get_gene_modules with non-existent gene."""
        result = get_gene_modules("nonexistent_gene_123")

        assert "error" in result
        assert "No modules found" in result["error"]

    def test_get_module_genes(self, loaded_metadata_registry):
        """Test get_module_genes function."""
        # Data is already loaded via fixture
        # Test with a known module if modules exist
        if loaded_metadata_registry.module_to_genes:
            module_id = list(loaded_metadata_registry.module_to_genes.keys())[0]
            result = get_module_genes(module_id)

            if "error" not in result:
                assert "data" in result
                data = result["data"]
                assert "module" in data
                assert "genes" in data
                assert "gene_count" in data

    def test_search_modules(self, loaded_metadata_registry):
        """Test search_modules function."""
        # Data is already loaded via fixture
        # Use a common term likely to find results
        result = search_modules("transport", 5)

        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        assert "modules" in data
        if data["modules"]:
            module_entry = data["modules"][0]
            assert "module" in module_entry
            assert "genes" in module_entry
            assert "gene_count" in module_entry

    def test_get_all_modules(self, loaded_metadata_registry):
        """Test get_all_modules function."""
        # Data is already loaded via fixture
        result = get_all_modules()

        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        assert "modules" in data
        if data["modules"]:
            module = data["modules"][0]
            assert "module_id" in module
            assert "name" in module
            assert "category" in module
            assert "count" in module


class TestComplexAnalysisFunctions:
    """Test cases for complex analysis functions with realistic scenarios."""

    def test_find_essential_genes_with_filter(self, loaded_metadata_registry):
        """Test find_essential_genes with condition filter."""
        # Data is already loaded via fixture
        result = find_essential_genes(
            condition_filter="pH", min_fitness_threshold=0.5, limit=3
        )

        # Check standardized response structure
        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        assert "essential_genes" in data
        # Could be empty if no essential genes found for pH conditions
        if data["essential_genes"]:
            for gene_entry in data["essential_genes"]:
                assert "essential_in_conditions" in gene_entry
                # Check that conditions contain pH filter
                for condition in gene_entry["essential_in_conditions"]:
                    condition_name = condition["condition"]
                    assert "pH" in condition_name or "ph" in condition_name.lower()

    def test_find_growth_inhibitor_genes_with_filter(self, loaded_metadata_registry):
        """Test find_growth_inhibitor_genes with condition filter."""
        # Data is already loaded via fixture
        result = find_growth_inhibitor_genes(
            condition_filter="stress", max_fitness_threshold=-0.5, limit=3
        )

        # Check standardized response structure
        assert isinstance(result, dict)
        assert "data" in result
        data = result["data"]
        assert "inhibitory_genes" in data
        # Could be empty if no inhibitor genes found for stress conditions
        if data["inhibitory_genes"]:
            for gene_entry in data["inhibitory_genes"]:
                assert "inhibits_growth_in_conditions" in gene_entry
                # Check that conditions contain stress filter
                for condition in gene_entry["inhibits_growth_in_conditions"]:
                    condition_name = condition["condition"]
                    assert "stress" in condition_name.lower()

    def test_gene_fitness_categorization(self, loaded_metadata_registry):
        """Test that fitness values are correctly categorized."""
        # Data is already loaded via fixture
        # Test with a gene that has diverse fitness values
        if loaded_metadata_registry.genes:
            gene_id = list(loaded_metadata_registry.genes.keys())[0]
            result = analyze_gene_fitness(gene_id)

            if "error" not in result:
                analysis = result["analysis"]

                # Check that categorization logic is working
                inhibitory = analysis["conditions_where_gene_inhibits_growth"]
                essential = analysis["conditions_where_gene_is_essential"]
                neutral = analysis["neutral_conditions"]

                # All should be lists
                assert isinstance(inhibitory, list)
                assert isinstance(essential, list)
                assert isinstance(neutral, list)

                # Check fitness value ranges if conditions exist
                for condition in inhibitory:
                    assert condition["fitness"] > 0.5  # Positive = gene inhibits growth

                for condition in essential:
                    assert condition["fitness"] < -0.5  # Negative = gene is essential

                for condition in neutral:
                    assert -0.5 <= condition["fitness"] <= 0.5
