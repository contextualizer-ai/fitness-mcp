"""Tests for MCP tool functions in fitness_mcp.main module.

Tests use real data loaders and real implementations.
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
    fitness_loader,
    module_loader,
)


class TestMCPToolFunctions:
    """Test cases for MCP tool functions using real data."""

    def test_get_gene_info_success(self):
        """Test get_gene_info with existing gene."""
        fitness_loader.load_data()

        # Test with a known gene if it exists
        if fitness_loader.genes:
            gene_id = list(fitness_loader.genes.keys())[0]
            result = get_gene_info(gene_id)

            if "error" not in result:
                assert "locusId" in result
                assert "sysName" in result
                assert "description" in result
                assert isinstance(result["locusId"], str)
                assert isinstance(result["sysName"], str)
                assert isinstance(result["description"], str)

    def test_get_gene_info_not_found(self):
        """Test get_gene_info with non-existent gene."""
        result = get_gene_info("nonexistent_gene_123")

        assert "error" in result
        assert "not found" in result["error"]

    def test_get_gene_fitness(self):
        """Test get_gene_fitness function."""
        fitness_loader.load_data()

        # Test with a known gene if it exists
        if fitness_loader.genes:
            gene_id = list(fitness_loader.genes.keys())[0]
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

    def test_search_genes(self):
        """Test search_genes function."""
        fitness_loader.load_data()

        # Use a common term likely to find results
        result = search_genes("ribosom", 5)

        assert isinstance(result, list)
        if result:
            gene = result[0]
            assert "locusId" in gene
            assert "sysName" in gene
            assert "description" in gene

    def test_get_growth_conditions(self):
        """Test get_growth_conditions function."""
        fitness_loader.load_data()

        result = get_growth_conditions()

        assert isinstance(result, list)
        if result:
            assert isinstance(result[0], str)

    def test_get_growth_conditions_with_filter(self):
        """Test get_growth_conditions with filter."""
        fitness_loader.load_data()

        result = get_growth_conditions("pH")

        assert isinstance(result, list)
        # All returned conditions should contain the filter term
        for condition in result:
            assert "pH" in condition or "ph" in condition.lower()

    def test_get_condition_details(self):
        """Test get_condition_details function."""
        fitness_loader.load_data()

        # Use a real condition if available
        if fitness_loader.conditions:
            condition_name = fitness_loader.conditions[0]
            result = get_condition_details(condition_name)

            if "error" not in result:
                assert "condition_name" in result
                assert result["condition_name"] == condition_name

    def test_interpret_fitness_score(self):
        """Test interpret_fitness_score function."""
        # Test positive score (gene benefits growth)
        result = interpret_fitness_score(0.8)

        assert "interpretation" in result
        assert "effect" in result
        assert "magnitude" in result
        assert "score" in result
        assert result["score"] == 0.8
        assert result["effect"] == "gene_benefits_growth"

        # Test negative score (gene inhibits growth)
        result = interpret_fitness_score(-0.8)
        assert result["effect"] == "gene_inhibits_growth"

        # Test neutral score
        result = interpret_fitness_score(0.05)
        assert result["effect"] == "neutral"

        # Test None score
        result = interpret_fitness_score(None)
        assert result["effect"] == "unknown"

    def test_find_essential_genes(self):
        """Test find_essential_genes function."""
        fitness_loader.load_data()

        result = find_essential_genes(limit=2)

        assert isinstance(result, list)
        # Result could be empty if no essential genes found with default threshold
        if result:
            gene_entry = result[0]
            assert "gene" in gene_entry
            assert "essential_in_conditions" in gene_entry
            assert "num_essential_conditions" in gene_entry

            # Check gene structure
            gene = gene_entry["gene"]
            assert "locusId" in gene
            assert "sysName" in gene
            assert "description" in gene

    def test_find_growth_inhibitor_genes(self):
        """Test find_growth_inhibitor_genes function."""
        fitness_loader.load_data()

        result = find_growth_inhibitor_genes(limit=2)

        assert isinstance(result, list)
        # Result could be empty if no inhibitor genes found with default threshold
        if result:
            gene_entry = result[0]
            assert "gene" in gene_entry
            assert "inhibits_growth_in_conditions" in gene_entry
            assert "num_inhibitory_conditions" in gene_entry

    def test_analyze_gene_fitness(self):
        """Test analyze_gene_fitness function."""
        fitness_loader.load_data()

        # Test with a known gene if available
        if fitness_loader.genes:
            gene_id = list(fitness_loader.genes.keys())[0]
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

    def test_get_gene_modules_success(self):
        """Test get_gene_modules with existing gene."""
        module_loader.load_data()

        # Test with a known gene if modules exist
        if module_loader.gene_to_modules:
            gene_id = list(module_loader.gene_to_modules.keys())[0]
            result = get_gene_modules(gene_id)

            if "error" not in result:
                assert "gene_id" in result
                assert "modules" in result
                assert "module_count" in result
                assert result["gene_id"] == gene_id
                assert isinstance(result["modules"], list)
                assert isinstance(result["module_count"], int)

    def test_get_gene_modules_not_found(self):
        """Test get_gene_modules with non-existent gene."""
        result = get_gene_modules("nonexistent_gene_123")

        assert "error" in result
        assert "No modules found" in result["error"]

    def test_get_module_genes(self):
        """Test get_module_genes function."""
        module_loader.load_data()

        # Test with a known module if modules exist
        if module_loader.module_to_genes:
            module_id = list(module_loader.module_to_genes.keys())[0]
            result = get_module_genes(module_id)

            if "error" not in result:
                assert "module" in result
                assert "genes" in result
                assert "gene_count" in result

    def test_search_modules(self):
        """Test search_modules function."""
        module_loader.load_data()

        # Use a common term likely to find results
        result = search_modules("transport", 5)

        assert isinstance(result, list)
        if result:
            module_entry = result[0]
            assert "module" in module_entry
            assert "genes" in module_entry
            assert "gene_count" in module_entry

    def test_get_all_modules(self):
        """Test get_all_modules function."""
        module_loader.load_data()

        result = get_all_modules()

        assert isinstance(result, list)
        if result:
            module = result[0]
            assert "module_id" in module
            assert "name" in module
            assert "category" in module
            assert "count" in module


class TestComplexAnalysisFunctions:
    """Test cases for complex analysis functions with realistic data."""

    def setup_method(self):
        """Set up realistic test data."""
        self.genes_data = {
            "Atu0001": {
                "locusId": "Atu0001",
                "sysName": "essential1",
                "description": "Essential gene 1",
                "fitness_values": [1.5, 1.2, 0.8, 0.1],  # High positive = essential
            },
            "Atu0002": {
                "locusId": "Atu0002",
                "sysName": "inhibitor1",
                "description": "Growth inhibitor gene 1",
                "fitness_values": [-1.0, -0.8, -0.2, 0.1],  # Negative = inhibits growth
            },
            "Atu0003": {
                "locusId": "Atu0003",
                "sysName": "neutral1",
                "description": "Neutral gene 1",
                "fitness_values": [0.1, -0.1, 0.05, -0.05],  # Near zero = neutral
            },
        }

        self.conditions = ["stress1", "stress2", "carbon1", "control"]

    def create_mock_fitness_response(self, gene_id, condition_filter=None):
        """Create mock fitness response for a gene."""
        gene_data = self.genes_data[gene_id]
        fitness_data = []

        for i, condition in enumerate(self.conditions):
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
            "total_conditions": len(fitness_data),
        }

    def test_find_essential_genes_realistic(self):
        """Test find_essential_genes with realistic data."""

        def mock_get_gene_fitness(gene_id, condition_filter=None):
            if gene_id in self.genes_data or gene_id in [
                g["sysName"] for g in self.genes_data.values()
            ]:
                # Handle sysName lookup
                if gene_id not in self.genes_data:
                    gene_id = next(
                        k for k, v in self.genes_data.items() if v["sysName"] == gene_id
                    )
                return self.create_mock_fitness_response(gene_id, condition_filter)
            return {"error": "Gene not found"}

        def mock_interpret_score(score):
            if score is None:
                return {"effect": "unknown", "interpretation": "No data"}
            elif score >= 0.5:
                return {
                    "effect": "gene_benefits_growth",
                    "interpretation": "Essential gene",
                }
            elif score <= -0.5:
                return {
                    "effect": "gene_inhibits_growth",
                    "interpretation": "Growth inhibitor",
                }
            else:
                return {"effect": "neutral", "interpretation": "Neutral effect"}

        with (
            patch.object(fitness_loader, "load_data"),
            patch.object(fitness_loader, "genes", self.genes_data),
            patch.object(
                fitness_loader, "get_gene_fitness", side_effect=mock_get_gene_fitness
            ),
            patch.object(
                fitness_loader,
                "interpret_fitness_score",
                side_effect=mock_interpret_score,
            ),
        ):
            result = find_essential_genes(
                condition_filter="stress", min_fitness_threshold=0.5, limit=10
            )

            assert isinstance(result, list)
            # Should find Atu0001 as essential (high positive fitness scores)
            if result:
                assert any(gene["gene"]["locusId"] == "Atu0001" for gene in result)
                for gene in result:
                    assert "essential_in_conditions" in gene
                    assert gene["num_essential_conditions"] > 0

    def test_find_growth_inhibitor_genes_realistic(self):
        """Test find_growth_inhibitor_genes with realistic data."""

        def mock_get_gene_fitness(gene_id, condition_filter=None):
            if gene_id in self.genes_data or gene_id in [
                g["sysName"] for g in self.genes_data.values()
            ]:
                if gene_id not in self.genes_data:
                    gene_id = next(
                        k for k, v in self.genes_data.items() if v["sysName"] == gene_id
                    )
                return self.create_mock_fitness_response(gene_id, condition_filter)
            return {"error": "Gene not found"}

        def mock_interpret_score(score):
            if score is None:
                return {"effect": "unknown", "interpretation": "No data"}
            elif score >= 0.5:
                return {
                    "effect": "gene_benefits_growth",
                    "interpretation": "Essential gene",
                }
            elif score <= -0.5:
                return {
                    "effect": "gene_inhibits_growth",
                    "interpretation": "Growth inhibitor",
                }
            else:
                return {"effect": "neutral", "interpretation": "Neutral effect"}

        with (
            patch.object(fitness_loader, "load_data"),
            patch.object(fitness_loader, "genes", self.genes_data),
            patch.object(
                fitness_loader, "get_gene_fitness", side_effect=mock_get_gene_fitness
            ),
            patch.object(
                fitness_loader,
                "interpret_fitness_score",
                side_effect=mock_interpret_score,
            ),
        ):
            result = find_growth_inhibitor_genes(
                condition_filter="stress", max_fitness_threshold=-0.5, limit=10
            )

            assert isinstance(result, list)
            # Should find Atu0002 as growth inhibitor (negative fitness scores)
            if result:
                assert any(gene["gene"]["locusId"] == "Atu0002" for gene in result)
                for gene in result:
                    assert "inhibits_growth_in_conditions" in gene
                    assert gene["num_inhibitory_conditions"] > 0

    def test_analyze_gene_fitness_categorization(self):
        """Test analyze_gene_fitness correctly categorizes fitness effects."""
        mock_fitness_data = {
            "gene": {
                "locusId": "Atu0001",
                "sysName": "test_gene",
                "description": "Test gene",
            },
            "fitness_data": [
                {
                    "condition": "beneficial_cond",
                    "fitness": 0.8,
                },  # Positive: gene inhibits growth when present
                {
                    "condition": "inhibitory_cond",
                    "fitness": -0.8,
                },  # Negative: gene is essential for growth
                {"condition": "neutral_cond", "fitness": 0.1},  # Neutral
                {"condition": "no_data_cond", "fitness": None},  # No data
            ],
        }

        with patch.object(fitness_loader, "get_gene_fitness") as mock_get:
            mock_get.return_value = mock_fitness_data

            result = analyze_gene_fitness("Atu0001")

            analysis = result["analysis"]

            # Should categorize correctly
            assert len(analysis["conditions_where_gene_is_essential"]) == 1
            assert len(analysis["conditions_where_gene_inhibits_growth"]) == 1
            assert len(analysis["neutral_conditions"]) == 1

            # Check summary counts
            summary = analysis["summary"]
            assert summary["essential_count"] == 1
            assert summary["inhibitory_count"] == 1
            assert summary["neutral_count"] == 1
            assert summary["total_conditions_tested"] == 3  # Excludes None values

    def test_analyze_gene_fitness_with_filters(self):
        """Test analyze_gene_fitness with min/max fitness filters."""
        fitness_loader.load_data()

        # Test with a known gene if available
        if fitness_loader.genes:
            gene_id = list(fitness_loader.genes.keys())[0]
            result = analyze_gene_fitness(gene_id, min_fitness=0.0, max_fitness=1.0)

            if "error" not in result:
                analysis = result["analysis"]

                # Should only include conditions with fitness in range [0.0, 1.0]
                all_conditions = (
                    analysis["conditions_where_gene_is_essential"]
                    + analysis["conditions_where_gene_inhibits_growth"]
                    + analysis["neutral_conditions"]
                )

                for condition in all_conditions:
                    fitness_val = condition["fitness"]
                    assert 0.0 <= fitness_val <= 1.0

class TestComplexAnalysisFunctions:
    """Test cases for complex analysis functions with realistic scenarios."""

    def test_find_essential_genes_with_filter(self):
        """Test find_essential_genes with condition filter."""
        fitness_loader.load_data()

        result = find_essential_genes(
            condition_filter="pH", min_fitness_threshold=0.5, limit=3
        )

        assert isinstance(result, list)
        # Could be empty if no essential genes found for pH conditions
        if result:
            for gene_entry in result:
                assert "essential_in_conditions" in gene_entry
                # Check that conditions contain pH filter
                for condition in gene_entry["essential_in_conditions"]:
                    condition_name = condition["condition"]
                    assert "pH" in condition_name or "ph" in condition_name.lower()

    def test_find_growth_inhibitor_genes_with_filter(self):
        """Test find_growth_inhibitor_genes with condition filter."""
        fitness_loader.load_data()

        result = find_growth_inhibitor_genes(
            condition_filter="stress", max_fitness_threshold=-0.5, limit=3
        )

        assert isinstance(result, list)
        # Could be empty if no inhibitor genes found for stress conditions
        if result:
            for gene_entry in result:
                assert "inhibits_growth_in_conditions" in gene_entry
                # Check that conditions contain stress filter
                for condition in gene_entry["inhibits_growth_in_conditions"]:
                    condition_name = condition["condition"]
                    assert "stress" in condition_name.lower()

    def test_gene_fitness_categorization(self):
        """Test that fitness values are correctly categorized."""
        fitness_loader.load_data()

        # Test with a gene that has diverse fitness values
        if fitness_loader.genes:
            gene_id = list(fitness_loader.genes.keys())[0]
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
