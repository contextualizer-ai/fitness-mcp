################################################################################
# fitness_mcp/main.py
# This module provides a FastMCP wrapper for gene fitness data analysis
################################################################################
import csv
import os
from typing import Any, Dict, List, Optional, Union
from fastmcp import FastMCP


# DATA LOADING SECTION
class ModuleDataLoader:
    def __init__(self, modules_file: str = "data/RbTnSeq_modules_t1e-7.csv", 
                 meta_file: str = "data/module_meta.tsv"):
        """Initialize the module data loader.
        
        Args:
            modules_file: Path to the modules CSV file
            meta_file: Path to the module metadata TSV file
        """
        self.modules_file = modules_file
        self.meta_file = meta_file
        self.gene_to_modules = {}  # gene_id -> [module_info, ...]
        self.module_to_genes = {}  # module_id -> [gene_info, ...]
        self.module_meta = {}      # module_id -> {name, category, count}
        self.loaded = False
        
    def load_data(self):
        """Load the module data from files."""
        if self.loaded:
            return
            
        # Load module metadata first
        if not os.path.exists(self.meta_file):
            raise FileNotFoundError(f"Module metadata file not found: {self.meta_file}")
            
        with open(self.meta_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                module_id = int(row['module'])
                self.module_meta[module_id] = {
                    'module_id': module_id,
                    'name': row['name'],
                    'category': row['category'],
                    'count': int(row['count'])
                }
        
        # Load gene-module assignments
        if not os.path.exists(self.modules_file):
            raise FileNotFoundError(f"Modules file not found: {self.modules_file}")
            
        with open(self.modules_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                module_id = int(row['module'])
                locus_tag = row['locus_tag']
                gene_weight = float(row['gene_weight']) if row['gene_weight'] else 0.0
                
                gene_info = {
                    'locus_tag': locus_tag,
                    'module_id': module_id,
                    'gene_weight': gene_weight,
                    'product': row.get('product', ''),
                    'description': row.get('Description', ''),
                    'preferred_name': row.get('Preferred_name', ''),
                    'protein_names': row.get('Protein names', ''),
                    'go_terms': row.get('Gene Ontology (GO)', '')
                }
                
                # Add module metadata to gene info
                if module_id in self.module_meta:
                    gene_info.update({
                        'module_name': self.module_meta[module_id]['name'],
                        'module_category': self.module_meta[module_id]['category']
                    })
                
                # Index by gene
                if locus_tag not in self.gene_to_modules:
                    self.gene_to_modules[locus_tag] = []
                self.gene_to_modules[locus_tag].append(gene_info)
                
                # Index by module
                if module_id not in self.module_to_genes:
                    self.module_to_genes[module_id] = []
                self.module_to_genes[module_id].append(gene_info)
                
        self.loaded = True
        
    def get_modules_for_gene(self, gene_id: str) -> List[Dict[str, Any]]:
        """Get module information for a specific gene.
        
        Args:
            gene_id: Gene locus tag
            
        Returns:
            List of modules containing this gene
        """
        self.load_data()
        return self.gene_to_modules.get(gene_id, [])
        
    def get_genes_in_module(self, module_id: int) -> Dict[str, Any]:
        """Get all genes in a specific module.
        
        Args:
            module_id: Module ID number
            
        Returns:
            Dictionary with module info and gene list
        """
        self.load_data()
        
        if module_id not in self.module_meta:
            return {'error': f'Module {module_id} not found'}
            
        genes = self.module_to_genes.get(module_id, [])
        
        return {
            'module': self.module_meta[module_id],
            'genes': genes,
            'gene_count': len(genes)
        }
        
    def search_modules_by_name(self, query: str, limit: int = 10) -> List[Dict[str, Any]]:
        """Search modules by name or category.
        
        Args:
            query: Search term
            limit: Maximum results
            
        Returns:
            List of matching modules with gene lists
        """
        self.load_data()
        
        query = query.lower()
        matches = []
        
        for module_id, meta in self.module_meta.items():
            if (query in meta['name'].lower() or 
                query in meta['category'].lower()):
                
                genes = self.module_to_genes.get(module_id, [])
                matches.append({
                    'module': meta,
                    'genes': genes,
                    'gene_count': len(genes)
                })
                
                if len(matches) >= limit:
                    break
                    
        return matches
        
    def get_all_modules(self) -> List[Dict[str, Any]]:
        """Get list of all modules with basic info.
        
        Returns:
            List of all modules
        """
        self.load_data()
        return list(self.module_meta.values())


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


# Global data loader instances
fitness_loader = FitnessDataLoader()
module_loader = ModuleDataLoader()


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


def get_gene_modules(
    gene_id: str
) -> Dict[str, Any]:
    """
    Get module information for a specific gene/locus.

    Args:
        gene_id: Gene locus tag (e.g., 'Atu0001')

    Returns:
        Dict containing modules that include this gene
    """
    modules = module_loader.get_modules_for_gene(gene_id)
    
    if not modules:
        return {'error': f'No modules found for gene {gene_id}'}
    
    return {
        'gene_id': gene_id,
        'modules': modules,
        'module_count': len(modules)
    }


def get_module_genes(
    module_id: int
) -> Dict[str, Any]:
    """
    Get all genes in a specific module.

    Args:
        module_id: Module ID number

    Returns:
        Dict containing module info and all genes in the module
    """
    return module_loader.get_genes_in_module(module_id)


def search_modules(
    query: str,
    limit: int = 10
) -> List[Dict[str, Any]]:
    """
    Search for modules by name or category.

    Args:
        query: Search term to match against module names or categories
        limit: Maximum number of results to return (default: 10)

    Returns:
        List of matching modules with their genes
    """
    return module_loader.search_modules_by_name(query, limit)


def get_all_modules() -> List[Dict[str, Any]]:
    """
    Get list of all available modules.

    Returns:
        List of all modules with basic information
    """
    return module_loader.get_all_modules()


# MAIN SECTION
# Create the FastMCP instance
mcp = FastMCP("fitness_mcp")

# Register all tools
mcp.tool(get_gene_info)
mcp.tool(get_gene_fitness)
mcp.tool(search_genes)
mcp.tool(get_growth_conditions)
mcp.tool(analyze_gene_fitness)
mcp.tool(get_gene_modules)
mcp.tool(get_module_genes)
mcp.tool(search_modules)
mcp.tool(get_all_modules)


def main():
    """Main entry point for the application."""
    mcp.run()


if __name__ == "__main__":
    main()
