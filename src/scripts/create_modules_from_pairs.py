import pandas as pd  # type: ignore
from collections import defaultdict


def create_modules_from_fitness_pairs(
    file_path, target_gene, threshold=2.0
):
    # Read the fitness pairs data
    df = pd.read_csv(
        file_path, 
        sep='\t', 
        names=['gene_id', 'condition_id', 'value'],
        dtype={'value': float},
        low_memory=False,
        skiprows=1  # Skip the header row
    )
    
    # Convert value to numeric, coercing errors to NaN
    df['value'] = pd.to_numeric(df['value'], errors='coerce')
    
    # Filter for significant fitness values (absolute value > threshold)
    significant_pairs = df[
        (df['gene_id'] == target_gene) & 
        (df['value'].notna()) & 
        (abs(df['value']) > threshold)
    ]
    
    # Group by conditions and find co-occurring genes
    condition_gene_map = defaultdict(set)
    for _, row in df[
        (df['condition_id'].isin(significant_pairs['condition_id'])) & 
        (abs(df['value']) > threshold)
    ].iterrows():
        condition_gene_map[row['condition_id']].add(row['gene_id'])
    
    # Create modules based on shared conditions
    modules = []
    processed_genes = set()
    
    for condition, genes in condition_gene_map.items():
        if target_gene not in genes:
            continue
        
        # Create a module for this condition
        module = set(genes)
        modules.append(module)
        processed_genes.update(module)
    
    # Merge overlapping modules
    merged_modules = []
    while modules:
        current_module = modules.pop(0)
        
        # Check for overlaps
        overlapping = [
            i for i, other_module in enumerate(modules) 
            if len(current_module.intersection(other_module)) > 0
        ]
        
        # Merge overlapping modules
        if overlapping:
            for idx in reversed(overlapping):
                current_module.update(modules.pop(idx))
        
        merged_modules.append(current_module)
    
    # Filter and format modules
    final_modules = []
    for module in merged_modules:
        if target_gene in module:
            final_modules.append(sorted(list(module)))
    
    return final_modules


# Run the analysis
file_path = 'data/fit_t_pairs_threshold_2_long.tab'
target_gene = 'Atu3150'
modules = create_modules_from_fitness_pairs(file_path, target_gene)

# Print results
print(f"Modules containing {target_gene}:")
for i, module in enumerate(modules, 1):
    print(f"Module {i}: {module}")
    print(f"Module size: {len(module)} genes")
    print()
