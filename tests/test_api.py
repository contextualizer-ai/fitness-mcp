"""Basic API tests for fitness_mcp."""

import pytest
from unittest.mock import patch, MagicMock

from src.fitness_mcp.main import (
    get_gene_info,
    get_gene_fitness,
    search_genes,
    get_growth_conditions,
    interpret_fitness_score,
    fitness_loader
)


def test_reality():
    """Basic sanity check."""
    assert 1 == 1


def test_gene_info_api_basic():
    """Test basic gene info API functionality."""
    with patch.object(fitness_loader, 'get_gene_info') as mock_get:
        mock_get.return_value = {
            "locusId": "Atu0001",
            "sysName": "testGene",
            "description": "Test gene description"
        }
        
        result = get_gene_info("Atu0001")
        
        assert "locusId" in result
        assert "sysName" in result
        assert "description" in result
        assert result["locusId"] == "Atu0001"


def test_gene_fitness_api_basic():
    """Test basic gene fitness API functionality."""
    mock_response = {
        "gene": {
            "locusId": "Atu0001",
            "sysName": "testGene", 
            "description": "Test gene"
        },
        "fitness_data": [
            {"condition": "test_cond", "fitness": 0.5}
        ],
        "total_conditions": 1
    }
    
    with patch.object(fitness_loader, 'get_gene_fitness') as mock_get:
        mock_get.return_value = mock_response
        
        result = get_gene_fitness("Atu0001")
        
        assert "gene" in result
        assert "fitness_data" in result
        assert len(result["fitness_data"]) == 1


def test_search_genes_api_basic():
    """Test basic search genes API functionality."""
    mock_results = [
        {
            "locusId": "Atu0001",
            "sysName": "gene1",
            "description": "Test gene 1"
        },
        {
            "locusId": "Atu0002", 
            "sysName": "gene2",
            "description": "Test gene 2"
        }
    ]
    
    with patch.object(fitness_loader, 'search_genes') as mock_search:
        mock_search.return_value = mock_results
        
        result = search_genes("test")
        
        assert isinstance(result, list)
        assert len(result) == 2
        assert all("locusId" in gene for gene in result)


def test_growth_conditions_api_basic():
    """Test basic growth conditions API functionality."""
    mock_conditions = ["condition1", "condition2", "condition3"]
    
    with patch.object(fitness_loader, 'get_conditions') as mock_get:
        mock_get.return_value = mock_conditions
        
        result = get_growth_conditions()
        
        assert isinstance(result, list)
        assert len(result) == 3


def test_interpret_fitness_score_api_basic():
    """Test basic fitness score interpretation API."""
    mock_interpretation = {
        "interpretation": "Gene knockout reduces fitness",
        "effect": "gene_benefits_growth", 
        "magnitude": "moderate",
        "score": 0.5
    }
    
    with patch.object(fitness_loader, 'interpret_fitness_score') as mock_interpret:
        mock_interpret.return_value = mock_interpretation
        
        result = interpret_fitness_score(0.5)
        
        assert "interpretation" in result
        assert "effect" in result
        assert "magnitude" in result
        assert result["score"] == 0.5
