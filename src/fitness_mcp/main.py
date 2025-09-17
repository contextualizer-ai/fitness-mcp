################################################################################
# fitness_mcp/main.py
# This module provides a FastMCP wrapper for gene fitness data analysis
################################################################################
import csv
import os
from typing import Any, Dict, List, Optional, Union
from fastmcp import FastMCP


# DATA LOADING SECTION
class FitnessDataLoader:
    def __init__(self, data_file: str = "data/fit_t.tab"):
        """Initialize the fitness data loader.
        
        Args:
            data_file: Path to the fitness data file
        """
        self.data_file = data_file
        self.genes = {}
        self.conditions = []
        self.loaded = False
        
    def load_data(self):
        """Load the fitness data from the tab-separated file."""
        if self.loaded:
            return
            
        if not os.path.exists(self.data_file):
            raise FileNotFoundError(f"Fitness data file not found: {self.data_file}")
            
        with open(self.data_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            
            # Read header (condition names)
            header = next(reader)
            self.conditions = header[3:]  # Skip locusId, sysName, desc columns
            
            # Read gene data
            for row in reader:
                if len(row) < 4:  # Skip incomplete rows
                    continue
                    
                locus_id = row[0]
                sys_name = row[1] 
                description = row[2]
                fitness_values = [float(val) if val and val != 'NA' else None 
                                for val in row[3:]]
                
                self.genes[locus_id] = {
                    'locusId': locus_id,
                    'sysName': sys_name,
                    'description': description,
                    'fitness_values': fitness_values
                }
                
                # Also index by sysName if different from locusId
                if sys_name and sys_name != locus_id:
                    self.genes[sys_name] = self.genes[locus_id]
                    
        self.loaded = True
        
    def get_gene_info(self, gene_id: str) -> Optional[Dict[str, Any]]:
        """Get basic information about a gene.
        
        Args:
            gene_id: Gene locus ID or system name
            
        Returns:
            Dictionary with gene information or None if not found
        """
        self.load_data()
        return self.genes.get(gene_id)
        
    def get_gene_fitness(self, gene_id: str, condition_filter: Optional[str] = None) -> Dict[str, Any]:
        """Get fitness data for a specific gene across growth conditions.
        
        Args:
            gene_id: Gene locus ID or system name
            condition_filter: Optional filter to match condition names (case-insensitive)
            
        Returns:
            Dictionary with gene info and fitness data
        """
        self.load_data()
        
        gene_info = self.genes.get(gene_id)
        if not gene_info:
            return {'error': f'Gene {gene_id} not found'}
            
        # Filter conditions if requested
        if condition_filter:
            condition_filter = condition_filter.lower()
            matching_conditions = []
            matching_values = []
            
            for i, condition in enumerate(self.conditions):
                if condition_filter in condition.lower():
                    matching_conditions.append(condition)
                    if i < len(gene_info['fitness_values']):
                        matching_values.append(gene_info['fitness_values'][i])
                    else:
                        matching_values.append(None)
        else:
            matching_conditions = self.conditions
            matching_values = gene_info['fitness_values']
            
        # Create condition-value pairs
        fitness_data = []
        for condition, value in zip(matching_conditions, matching_values):
            fitness_data.append({
                'condition': condition,
                'fitness': value
            })
            
        return {
            'gene': {
                'locusId': gene_info['locusId'],
                'sysName': gene_info['sysName'], 
                'description': gene_info['description']
            },
            'fitness_data': fitness_data,
            'total_conditions': len(fitness_data)
        }
        
    def search_genes(self, query: str, limit: int = 10) -> List[Dict[str, Any]]:
        """Search for genes by name or description.
        
        Args:
            query: Search term to match against gene names or descriptions
            limit: Maximum number of results to return
            
        Returns:
            List of matching gene information
        """
        self.load_data()
        
        query = query.lower()
        matches = []
        seen_locus_ids = set()
        
        for gene_id, gene_info in self.genes.items():
            # Skip duplicates (since we index by both locusId and sysName)
            if gene_info['locusId'] in seen_locus_ids:
                continue
                
            # Check if query matches locusId, sysName, or description
            if (query in gene_info['locusId'].lower() or 
                query in gene_info['sysName'].lower() or
                query in gene_info['description'].lower()):
                
                matches.append({
                    'locusId': gene_info['locusId'],
                    'sysName': gene_info['sysName'],
                    'description': gene_info['description']
                })
                seen_locus_ids.add(gene_info['locusId'])
                
                if len(matches) >= limit:
                    break
                    
        return matches
        
    def get_conditions(self, condition_filter: Optional[str] = None) -> List[str]:
        """Get list of available growth conditions.
        
        Args:
            condition_filter: Optional filter to match condition names
            
        Returns:
            List of condition names
        """
        self.load_data()
        
        if condition_filter:
            condition_filter = condition_filter.lower()
            return [cond for cond in self.conditions if condition_filter in cond.lower()]
        else:
            return self.conditions


# Global data loader instance
fitness_loader = FitnessDataLoader()


# MCP TOOL SECTION
def get_gene_info(
    gene_id: str
) -> Dict[str, Any]:
    """
    Get basic information about a gene.

    Args:
        gene_id: Gene locus ID (e.g., 'Atu0001') or system name

    Returns:
        Dict containing gene information including locusId, sysName, and description
    """
    gene_info = fitness_loader.get_gene_info(gene_id)
    if gene_info:
        return {
            'locusId': gene_info['locusId'],
            'sysName': gene_info['sysName'],
            'description': gene_info['description']
        }
    else:
        return {'error': f'Gene {gene_id} not found'}


def get_gene_fitness(
    gene_id: str,
    condition_filter: Optional[str] = None
) -> Dict[str, Any]:
    """
    Get fitness data for a gene across different growth conditions.

    Args:
        gene_id: Gene locus ID (e.g., 'Atu0001') or system name
        condition_filter: Optional filter to match specific conditions (e.g., 'glucose', 'stress', 'antibiotic')

    Returns:
        Dict containing gene info and fitness data for matching conditions
    """
    return fitness_loader.get_gene_fitness(gene_id, condition_filter)


def search_genes(
    query: str,
    limit: int = 10
) -> List[Dict[str, Any]]:
    """
    Search for genes by name or description.

    Args:
        query: Search term to match against gene names or descriptions
        limit: Maximum number of results to return (default: 10)

    Returns:
        List of dictionaries containing matching gene information
    """
    return fitness_loader.search_genes(query, limit)


def get_growth_conditions(
    condition_filter: Optional[str] = None
) -> List[str]:
    """
    Get list of available growth conditions.

    Args:
        condition_filter: Optional filter to match condition names (e.g., 'LB', 'glucose', 'metal')

    Returns:
        List of condition names
    """
    return fitness_loader.get_conditions(condition_filter)


def analyze_gene_fitness(
    gene_id: str,
    min_fitness: Optional[float] = None,
    max_fitness: Optional[float] = None
) -> Dict[str, Any]:
    """
    Analyze fitness effects for a gene with optional filtering by fitness values.

    Args:
        gene_id: Gene locus ID or system name
        min_fitness: Minimum fitness value to include (for identifying beneficial conditions)
        max_fitness: Maximum fitness value to include (for identifying detrimental conditions)

    Returns:
        Dict with gene info and categorized fitness data
    """
    fitness_data = fitness_loader.get_gene_fitness(gene_id)
    
    if 'error' in fitness_data:
        return fitness_data
        
    # Categorize conditions by fitness effect
    beneficial = []  # High fitness (gene knockout is beneficial)
    detrimental = []  # Low fitness (gene knockout is detrimental/essential)
    neutral = []  # Near-zero fitness (no significant effect)
    
    for item in fitness_data['fitness_data']:
        fitness_val = item['fitness']
        if fitness_val is None:
            continue
            
        # Apply filters if specified
        if min_fitness is not None and fitness_val < min_fitness:
            continue
        if max_fitness is not None and fitness_val > max_fitness:
            continue
            
        if fitness_val > 0.5:  # Significantly beneficial
            beneficial.append(item)
        elif fitness_val < -0.5:  # Significantly detrimental
            detrimental.append(item)
        else:  # Neutral effect
            neutral.append(item)
            
    return {
        'gene': fitness_data['gene'],
        'analysis': {
            'beneficial_conditions': beneficial,
            'detrimental_conditions': detrimental, 
            'neutral_conditions': neutral,
            'summary': {
                'total_conditions_tested': len([x for x in fitness_data['fitness_data'] if x['fitness'] is not None]),
                'beneficial_count': len(beneficial),
                'detrimental_count': len(detrimental),
                'neutral_count': len(neutral)
            }
        }
    }


# MAIN SECTION
# Create the FastMCP instance
mcp = FastMCP("fitness_mcp")

# Register all tools
mcp.tool(get_gene_info)
mcp.tool(get_gene_fitness)
mcp.tool(search_genes)
mcp.tool(get_growth_conditions)
mcp.tool(analyze_gene_fitness)


def main():
    """Main entry point for the application."""
    mcp.run()


if __name__ == "__main__":
    main()
