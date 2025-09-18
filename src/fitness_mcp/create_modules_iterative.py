#!/usr/bin/env python3
"""
Create modules from pairs data using iterative approach.

Starting from a single gene-condition pair, alternately add:
1. A gene that shares the current condition(s)
2. A condition that is shared by all current genes

This creates tight, coherent modules rather than large submatrices.
"""

import csv
import random
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
from collections import defaultdict


class IterativeModuleBuilder:
    def __init__(self, pairs_file: str):
        """Initialize with pairs data file."""
        self.pairs_file = pairs_file
        self.gene_to_conditions: Dict[str, Set[str]] = defaultdict(set)
        self.condition_to_genes: Dict[str, Set[str]] = defaultdict(set)
        self.all_pairs: Set[Tuple[str, str]] = set()
        self.load_pairs_data()
    
    def load_pairs_data(self):
        """Load gene-condition pairs from the data file."""
        with open(self.pairs_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                gene_id = row['gene_id']
                condition_id = row['condition_id']
                
                self.gene_to_conditions[gene_id].add(condition_id)
                self.condition_to_genes[condition_id].add(gene_id)
                self.all_pairs.add((gene_id, condition_id))
    
    def get_genes_with_all_conditions(self, conditions: Set[str]) -> Set[str]:
        """Find genes that have pairs with ALL the given conditions."""
        if not conditions:
            return set()
        
        # Start with genes from the first condition
        result = self.condition_to_genes[next(iter(conditions))].copy()
        
        # Intersect with genes from each remaining condition
        for condition in conditions:
            result &= self.condition_to_genes[condition]
        
        return result
    
    def get_conditions_with_all_genes(self, genes: Set[str]) -> Set[str]:
        """Find conditions that have pairs with ALL the given genes."""
        if not genes:
            return set()
        
        # Start with conditions from the first gene
        result = self.gene_to_conditions[next(iter(genes))].copy()
        
        # Intersect with conditions from each remaining gene
        for gene in genes:
            result &= self.gene_to_conditions[gene]
        
        return result
    
    def build_module(self, start_gene: str, start_condition: str, max_size: int = 10) -> Dict:
        """
        Build a module starting from a gene-condition pair.
        
        Args:
            start_gene: Starting gene
            start_condition: Starting condition
            max_size: Maximum module size to prevent runaway growth
            
        Returns:
            Dict with module genes, conditions, and statistics
        """
        # Verify the starting pair exists
        if (start_gene, start_condition) not in self.all_pairs:
            return {
                "error": f"Starting pair ({start_gene}, {start_condition}) not found in data"
            }
        
        module_genes = {start_gene}
        module_conditions = {start_condition}
        
        iteration = 0
        max_iterations = max_size * 2  # Prevent infinite loops
        
        history = []
        
        while iteration < max_iterations:
            iteration += 1
            
            # Alternate between adding genes and conditions
            if iteration % 2 == 1:
                # Odd iteration: Add a gene that has pairs with ALL current conditions
                candidate_genes = self.get_genes_with_all_conditions(module_conditions)
                candidate_genes -= module_genes  # Remove genes already in module
                
                if not candidate_genes:
                    history.append(f"Iteration {iteration}: No new genes found with all conditions")
                    break
                
                # Randomly select a new gene
                new_gene = random.choice(list(candidate_genes))
                module_genes.add(new_gene)
                history.append(f"Iteration {iteration}: Added gene {new_gene}")
                
            else:
                # Even iteration: Add a condition that has pairs with ALL current genes
                candidate_conditions = self.get_conditions_with_all_genes(module_genes)
                candidate_conditions -= module_conditions  # Remove conditions already in module
                
                if not candidate_conditions:
                    history.append(f"Iteration {iteration}: No new conditions found with all genes")
                    break
                
                # Randomly select a new condition
                new_condition = random.choice(list(candidate_conditions))
                module_conditions.add(new_condition)
                history.append(f"Iteration {iteration}: Added condition {new_condition}")
            
            # Stop if we've reached max size - this indicates runaway growth
            if len(module_genes) + len(module_conditions) >= max_size:
                history.append(f"Iteration {iteration}: Reached maximum size limit - runaway growth detected")
                return {
                    "error": f"Module building terminated due to runaway growth (>{max_size} elements)",
                    "reason": "Large modules indicate overly broad relationships rather than tight functional clusters",
                    "start_gene": start_gene,
                    "start_condition": start_condition,
                    "partial_build": {
                        "genes_found": len(module_genes),
                        "conditions_found": len(module_conditions),
                        "iterations": iteration,
                        "history": history[-3:]  # Show last few steps
                    }
                }
        
        # Require at least 2 genes and 2 conditions for a meaningful module
        if len(module_genes) < 2 or len(module_conditions) < 2:
            return {
                "error": "Module too small to be meaningful",
                "reason": f"Found only {len(module_genes)} genes and {len(module_conditions)} conditions",
                "start_gene": start_gene,
                "start_condition": start_condition,
                "iterations": iteration,
                "history": history
            }
        
        # Calculate module density
        density = len(module_genes) * len(module_conditions) / (len(module_genes) + len(module_conditions))
        
        return {
            "start_gene": start_gene,
            "start_condition": start_condition,
            "module_genes": sorted(module_genes),
            "module_conditions": sorted(module_conditions),
            "num_genes": len(module_genes),
            "num_conditions": len(module_conditions),
            "total_pairs": len(module_genes) * len(module_conditions),
            "density": round(density, 2),
            "iterations": iteration,
            "stopped_naturally": True,
            "history": history
        }
    
    def build_multiple_modules(self, num_modules: int = 5, max_size: int = 20, 
                              seed: Optional[int] = None) -> List[Dict]:
        """
        Build multiple modules using random starting points.
        
        Args:
            num_modules: Number of modules to build
            max_size: Maximum size per module
            seed: Random seed for reproducibility
            
        Returns:
            List of module dictionaries
        """
        if seed is not None:
            random.seed(seed)
        
        modules = []
        used_starts: Set[Tuple[str, str]] = set()
        
        for i in range(num_modules):
            # Choose a random starting pair
            available_pairs = list(self.all_pairs - used_starts)
            if not available_pairs:
                break
            
            start_gene, start_condition = random.choice(available_pairs)
            used_starts.add((start_gene, start_condition))
            
            module = self.build_module(start_gene, start_condition, max_size)
            if "error" not in module:
                modules.append(module)
        
        return modules


def main():
    """Main function to demonstrate the iterative module building."""
    project_root = Path(__file__).parent.parent.parent
    pairs_file = project_root / 'data' / 'fit_t_pairs_threshold_2_long.tab'
    
    if not pairs_file.exists():
        print(f"Error: Pairs file not found at {pairs_file}")
        return
    
    print("🧬 Building modules using iterative approach...")
    print("=" * 50)
    
    builder = IterativeModuleBuilder(str(pairs_file))
    
    # Build a specific module starting from a known interesting gene
    print("Building module starting from Atu3150 (lactose transporter)...")
    
    # Find a condition where Atu3150 has a significant fitness value
    atu3150_conditions = builder.gene_to_conditions.get('Atu3150', set())
    if atu3150_conditions:
        start_condition = random.choice(list(atu3150_conditions))
        module = builder.build_module('Atu3150', start_condition, max_size=15)
        
        if "error" not in module:
            print(f"\nSuccessful module starting from Atu3150 + {start_condition}:")
            print(f"Genes ({module['num_genes']}): {module['module_genes']}")
            print(f"Conditions ({module['num_conditions']}): {module['module_conditions']}")
            print(f"Total possible pairs: {module['total_pairs']}")
            print(f"Density: {module['density']}")
            print(f"Iterations: {module['iterations']}")
            print("\nHistory:")
            for step in module['history']:
                print(f"  {step}")
        else:
            print(f"❌ Failed to build module: {module['error']}")
            print(f"   Reason: {module['reason']}")
            if 'partial_build' in module:
                print(f"   Partial progress: {module['partial_build']['genes_found']} genes, {module['partial_build']['conditions_found']} conditions")
    else:
        print("Atu3150 not found in pairs data")
    
    print("\n" + "=" * 50)
    print("Building 3 random modules...")
    
    # Build multiple random modules
    modules = builder.build_multiple_modules(num_modules=5, max_size=8, seed=42)
    
    successful_modules = [m for m in modules if "error" not in m]
    failed_modules = [m for m in modules if "error" in m]
    
    print(f"Successfully built {len(successful_modules)} modules:")
    for i, module in enumerate(successful_modules, 1):
        print(f"\n✅ Module {i}:")
        print(f"  Start: {module['start_gene']} + {module['start_condition']}")
        print(f"  Genes ({module['num_genes']}): {', '.join(module['module_genes'][:5])}{'...' if module['num_genes'] > 5 else ''}")
        print(f"  Conditions ({module['num_conditions']}): {', '.join(module['module_conditions'][:3])}{'...' if module['num_conditions'] > 3 else ''}")
        print(f"  Total pairs: {module['total_pairs']}")
        print(f"  Density: {module['density']}")
        print(f"  Iterations: {module['iterations']}")
    
    if failed_modules:
        print(f"\n❌ Failed to build {len(failed_modules)} modules (runaway growth or too small)")
        for i, module in enumerate(failed_modules, 1):
            print(f"  {i}. {module['start_gene']} + {module['start_condition']}: {module['error']}")


if __name__ == "__main__":
    main()