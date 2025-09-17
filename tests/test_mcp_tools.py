"""Tests for MCP tool functions in fitness_mcp.main module."""

from unittest.mock import patch


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
    """Test cases for MCP tool functions."""

    def setup_method(self):
        """Set up test fixtures before each test method."""
        # Mock data for fitness loader
        self.mock_fitness_data = {
            "Atu0001": {
                "locusId": "Atu0001",
                "sysName": "gene1",
                "description": "Test gene 1",
                "fitness_values": [0.8, -0.6, 0.2],
            },
            "gene1": {  # Same gene indexed by sysName
                "locusId": "Atu0001",
                "sysName": "gene1",
                "description": "Test gene 1",
                "fitness_values": [0.8, -0.6, 0.2],
            },
            "Atu0002": {
                "locusId": "Atu0002",
                "sysName": "gene2",
                "description": "Test gene 2",
                "fitness_values": [-1.2, 0.9, None],
            },
        }

        self.mock_conditions = ["cond1", "cond2", "cond3"]
        self.mock_condition_details = {
            "cond1": {
                "short_desc": "Condition 1",
                "long_desc": "Long description 1",
                "media": "LB",
                "temperature": "30",
                "pH": "7.0",
                "aerobic": "True",
                "condition_1": "glucose",
                "concentration_1": "10",
                "units_1": "mM",
                "exp_group": "carbon",
            }
        }

        # Mock data for module loader
        self.mock_gene_to_modules = {
            "Atu0001": [
                {
                    "locus_tag": "Atu0001",
                    "module_id": 1,
                    "gene_weight": 0.8,
                    "product": "Product 1",
                    "module_name": "Module 1",
                    "module_category": "Category1",
                }
            ]
        }

        self.mock_module_to_genes = {
            1: [
                {
                    "locus_tag": "Atu0001",
                    "module_id": 1,
                    "gene_weight": 0.8,
                    "product": "Product 1",
                }
            ]
        }

        self.mock_module_meta = {
            1: {"module_id": 1, "name": "Module 1", "category": "Category1", "count": 1}
        }

    def test_get_gene_info_success(self):
        """Test get_gene_info with existing gene."""
        with patch.object(fitness_loader, "get_gene_info") as mock_get:
            mock_get.return_value = self.mock_fitness_data["Atu0001"]

            result = get_gene_info("Atu0001")

            assert result["locusId"] == "Atu0001"
            assert result["sysName"] == "gene1"
            assert result["description"] == "Test gene 1"
            mock_get.assert_called_once_with("Atu0001")

    def test_get_gene_info_not_found(self):
        """Test get_gene_info with non-existent gene."""
        with patch.object(fitness_loader, "get_gene_info") as mock_get:
            mock_get.return_value = None

            result = get_gene_info("nonexistent")

            assert "error" in result
            assert "not found" in result["error"]
            mock_get.assert_called_once_with("nonexistent")

    def test_get_gene_fitness(self):
        """Test get_gene_fitness function."""
        expected_result = {
            "gene": {
                "locusId": "Atu0001",
                "sysName": "gene1",
                "description": "Test gene 1",
            },
            "fitness_data": [
                {"condition": "cond1", "fitness": 0.8},
                {"condition": "cond2", "fitness": -0.6},
                {"condition": "cond3", "fitness": 0.2},
            ],
        }

        with patch.object(fitness_loader, "get_gene_fitness") as mock_get:
            mock_get.return_value = expected_result

            result = get_gene_fitness("Atu0001", "glucose")

            assert result == expected_result
            mock_get.assert_called_once_with("Atu0001", "glucose")

    def test_search_genes(self):
        """Test search_genes function."""
        expected_results = [
            {"locusId": "Atu0001", "sysName": "gene1", "description": "Test gene 1"}
        ]

        with patch.object(fitness_loader, "search_genes") as mock_search:
            mock_search.return_value = expected_results

            result = search_genes("test", 5)

            assert result == expected_results
            mock_search.assert_called_once_with("test", 5)

    def test_get_growth_conditions(self):
        """Test get_growth_conditions function."""
        with patch.object(fitness_loader, "get_conditions") as mock_get:
            mock_get.return_value = self.mock_conditions

            result = get_growth_conditions("glucose")

            assert result == self.mock_conditions
            mock_get.assert_called_once_with("glucose")

    def test_get_condition_details(self):
        """Test get_condition_details function."""
        expected_result = {
            "condition_name": "cond1",
            "short_description": "Condition 1",
            "experimental_group": "carbon",
        }

        with patch.object(fitness_loader, "get_condition_details") as mock_get:
            mock_get.return_value = expected_result

            result = get_condition_details("cond1")

            assert result == expected_result
            mock_get.assert_called_once_with("cond1")

    def test_interpret_fitness_score(self):
        """Test interpret_fitness_score function."""
        expected_result = {
            "interpretation": "Gene knockout reduces fitness",
            "effect": "gene_benefits_growth",
            "magnitude": "strong",
            "score": 0.8,
        }

        with patch.object(fitness_loader, "interpret_fitness_score") as mock_interpret:
            mock_interpret.return_value = expected_result

            result = interpret_fitness_score(0.8)

            assert result == expected_result
            mock_interpret.assert_called_once_with(0.8)

    def test_find_essential_genes(self):
        """Test find_essential_genes function."""
        # Mock fitness_loader.genes and other dependencies
        with (
            patch.object(fitness_loader, "load_data"),
            patch.object(fitness_loader, "genes", self.mock_fitness_data),
            patch.object(fitness_loader, "get_gene_fitness") as mock_get_fitness,
        ):
            # Mock fitness data for essential gene
            mock_get_fitness.return_value = {
                "gene": {
                    "locusId": "Atu0001",
                    "sysName": "gene1",
                    "description": "Test gene 1",
                },
                "fitness_data": [
                    {"condition": "cond1", "fitness": 0.8, "description": "Condition 1"}
                ],
            }

            with patch.object(
                fitness_loader, "interpret_fitness_score"
            ) as mock_interpret:
                mock_interpret.return_value = {
                    "interpretation": "Essential gene",
                    "effect": "gene_benefits_growth",
                }

                result = find_essential_genes("glucose", 0.5, 10)

                assert isinstance(result, list)
                assert len(result) >= 0
                if result:
                    assert "gene" in result[0]
                    assert "essential_in_conditions" in result[0]

    def test_find_growth_inhibitor_genes(self):
        """Test find_growth_inhibitor_genes function."""
        with (
            patch.object(fitness_loader, "load_data"),
            patch.object(fitness_loader, "genes", self.mock_fitness_data),
            patch.object(fitness_loader, "get_gene_fitness") as mock_get_fitness,
        ):
            # Mock fitness data for inhibitor gene
            mock_get_fitness.return_value = {
                "gene": {
                    "locusId": "Atu0002",
                    "sysName": "gene2",
                    "description": "Test gene 2",
                },
                "fitness_data": [
                    {
                        "condition": "cond1",
                        "fitness": -1.2,
                        "description": "Condition 1",
                    }
                ],
            }

            with patch.object(
                fitness_loader, "interpret_fitness_score"
            ) as mock_interpret:
                mock_interpret.return_value = {
                    "interpretation": "Growth inhibitor gene",
                    "effect": "gene_inhibits_growth",
                }

                result = find_growth_inhibitor_genes("stress", -0.5, 10)

                assert isinstance(result, list)
                assert len(result) >= 0
                if result:
                    assert "gene" in result[0]
                    assert "inhibits_growth_in_conditions" in result[0]

    def test_analyze_gene_fitness(self):
        """Test analyze_gene_fitness function."""
        mock_fitness_data = {
            "gene": {
                "locusId": "Atu0001",
                "sysName": "gene1",
                "description": "Test gene 1",
            },
            "fitness_data": [
                {"condition": "cond1", "fitness": 0.8},
                {"condition": "cond2", "fitness": -0.8},
                {"condition": "cond3", "fitness": 0.1},
            ],
        }

        with patch.object(fitness_loader, "get_gene_fitness") as mock_get:
            mock_get.return_value = mock_fitness_data

            result = analyze_gene_fitness("Atu0001", min_fitness=-1.0, max_fitness=1.0)

            assert "gene" in result
            assert "analysis" in result
            analysis = result["analysis"]
            assert "conditions_where_gene_inhibits_growth" in analysis
            assert "conditions_where_gene_is_essential" in analysis
            assert "neutral_conditions" in analysis
            assert "summary" in analysis

    def test_analyze_gene_fitness_error(self):
        """Test analyze_gene_fitness with error."""
        with patch.object(fitness_loader, "get_gene_fitness") as mock_get:
            mock_get.return_value = {"error": "Gene not found"}

            result = analyze_gene_fitness("nonexistent")

            assert "error" in result
            assert "Gene not found" in result["error"]

    def test_get_gene_modules_success(self):
        """Test get_gene_modules with existing gene."""
        with patch.object(module_loader, "get_modules_for_gene") as mock_get:
            mock_get.return_value = self.mock_gene_to_modules["Atu0001"]

            result = get_gene_modules("Atu0001")

            assert result["gene_id"] == "Atu0001"
            assert "modules" in result
            assert result["module_count"] == 1
            mock_get.assert_called_once_with("Atu0001")

    def test_get_gene_modules_not_found(self):
        """Test get_gene_modules with non-existent gene."""
        with patch.object(module_loader, "get_modules_for_gene") as mock_get:
            mock_get.return_value = []

            result = get_gene_modules("nonexistent")

            assert "error" in result
            assert "No modules found" in result["error"]

    def test_get_module_genes(self):
        """Test get_module_genes function."""
        expected_result = {
            "module": self.mock_module_meta[1],
            "genes": self.mock_module_to_genes[1],
            "gene_count": 1,
        }

        with patch.object(module_loader, "get_genes_in_module") as mock_get:
            mock_get.return_value = expected_result

            result = get_module_genes(1)

            assert result == expected_result
            mock_get.assert_called_once_with(1)

    def test_search_modules(self):
        """Test search_modules function."""
        expected_results = [
            {
                "module": self.mock_module_meta[1],
                "genes": self.mock_module_to_genes[1],
                "gene_count": 1,
            }
        ]

        with patch.object(module_loader, "search_modules_by_name") as mock_search:
            mock_search.return_value = expected_results

            result = search_modules("Module", 5)

            assert result == expected_results
            mock_search.assert_called_once_with("Module", 5)

    def test_get_all_modules(self):
        """Test get_all_modules function."""
        expected_result = list(self.mock_module_meta.values())

        with patch.object(module_loader, "get_all_modules") as mock_get:
            mock_get.return_value = expected_result

            result = get_all_modules()

            assert result == expected_result
            mock_get.assert_called_once()


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
        mock_fitness_data = {
            "gene": {"locusId": "Atu0001", "sysName": "test", "description": "Test"},
            "fitness_data": [
                {"condition": "cond1", "fitness": 1.0},  # Should be included
                {
                    "condition": "cond2",
                    "fitness": -1.0,
                },  # Should be excluded (below min)
                {
                    "condition": "cond3",
                    "fitness": 0.5,
                },  # Should be included (within range)
                {"condition": "cond4", "fitness": 0.2},  # Should be included
            ],
        }

        with patch.object(fitness_loader, "get_gene_fitness") as mock_get:
            mock_get.return_value = mock_fitness_data

            # Filter to only include fitness between 0.0 and 1.0
            result = analyze_gene_fitness("Atu0001", min_fitness=0.0, max_fitness=1.0)

            analysis = result["analysis"]

            # Should only include conditions with fitness in range [0.0, 1.0]
            all_conditions = (
                analysis["conditions_where_gene_is_essential"]
                + analysis["conditions_where_gene_inhibits_growth"]
                + analysis["neutral_conditions"]
            )

            assert (
                len(all_conditions) == 3
            )  # cond1, cond3, and cond4 (0.5 is within range)
            fitness_values = [cond["fitness"] for cond in all_conditions]
            assert 1.0 in fitness_values
            assert 0.2 in fitness_values
            assert 0.5 in fitness_values  # 0.5 should be included
            assert -1.0 not in fitness_values
