# Fitness MCP Development Guide

This document establishes development patterns and guidelines for the fitness-mcp project. Follow these patterns consistently for all development work.

## Project Overview

The Fitness MCP is a FastMCP server for analyzing gene fitness data from barcoded Agrobacterium mutant libraries. It efficiently handles a 90MB TSV dataset and provides scientific analysis tools for competitive fitness experiments.

## Key Development Commands

### Essential Commands
- Run MCP server: `uv run fitness-mcp`
- Test fitness protocol: `make test-fitness-protocol`
- Test gene analysis: `make test-gene-analysis`
- Claude demo: `make demo-atu3150-function`
- Full development cycle: `make all`

### GitHub Issue Workflow
For each GitHub issue, create a linked branch that will auto-close the issue when merged:

```bash
# Create and switch to issue branch (auto-links to GitHub issue)
git checkout -b issue-21-metadata-registry

# Work on the issue, making commits
git add .
git commit -m "Implement MetadataRegistry class

Addresses issue #21 by creating centralized metadata access.
🤖 Generated with [Claude Code](https://claude.ai/code)"

# When ready, push and create PR
git push -u origin issue-21-metadata-registry
gh pr create --title "Implement Centralized Metadata Registry" --body "Closes #21"
```

**Branch Naming Convention**: `issue-{number}-{short-description}`
- `issue-21-metadata-registry`
- `issue-22-eliminate-redundancy`
- `issue-23-consistent-naming`
- `issue-24-posix-tools`

## Architecture Principles

### Data Loading Best Practices
The MCP follows TSV best practices for large file handling:

```python
# ✅ Good: Single shared data loader
fitness_loader = FitnessDataLoader()  # Global instance shared across tools

# ✅ Good: Thread-safe operations
with self._lock:
    if not self._needs_reload():
        return

# ✅ Good: File change detection
def _needs_reload(self) -> bool:
    try:
        data_mtime = os.path.getmtime(self.data_file)
        exp_mtime = os.path.getmtime(self.exp_desc_file)
    except FileNotFoundError:
        return False
    return (not self.loaded or 
            data_mtime != self._mtime or 
            exp_mtime != self._exp_mtime)

# ✅ Good: LRU caching with data versioning
@lru_cache(maxsize=256)
def _cached_search_genes(self, query: str, limit: int, data_version: float):
    # Expensive operations cached with data file timestamp as key
```

### Anti-Patterns to Avoid

```python
# ❌ Bad: Loading data in every tool call
def search_genes(query):
    df = pd.read_csv("data/fit_t.tab", sep="\t")  # DON'T DO THIS
    return df[df['description'].str.contains(query)]

# ❌ Bad: No thread safety
def load_data(self):
    if not self.loaded:  # Race condition risk
        # load data without locks

# ❌ Bad: No caching
def search_genes(self, query: str):
    # Expensive search operation without caching
    matches = []
    for gene_id, gene_info in self.genes.items():  # Repeated work
        # ... expensive string matching
```

## Biological Data Standards

### Fitness Score Interpretation
**Critical**: Maintain correct biological interpretation throughout codebase:

```python
# ✅ Good: Correct interpretation in docstrings and code
"""
Fitness scores represent gene knockout effects:
- Negative values: Gene knockout IMPROVES fitness (gene normally inhibits growth)
- Positive values: Gene knockout REDUCES fitness (gene is beneficial/essential)
- Values near 0: Gene knockout has MINIMAL effect on fitness
"""

# ✅ Good: Biologically meaningful variable names
gene_inhibits_growth = []     # for negative fitness scores
gene_benefits_growth = []     # for positive fitness scores
neutral = []                  # for scores near zero

# ❌ Bad: Confusing or backwards interpretation
beneficial = []   # Ambiguous - beneficial to whom?
detrimental = []  # Confusing - could mean either direction
```

### Scientific Accuracy Requirements
- All docstrings must correctly explain biological meaning
- Function names should reflect biological interpretation
- Error messages should use proper scientific terminology
- Comments should explain the experimental context

## MCP Tool Design Patterns

### Tool Function Structure
Follow this pattern for all MCP tools:

```python
def tool_name(
    required_param: str,
    optional_param: Optional[str] = None
) -> Dict[str, Any]:
    """
    Brief description of what the tool does.

    Args:
        required_param: Description with example (e.g., 'Atu0001')
        optional_param: Description with examples (e.g., 'glucose', 'pH')

    Returns:
        Dict with structured results. Explain the data structure and biological meaning.
    """
    return fitness_loader.method_name(required_param, optional_param)
```

### Error Handling
```python
# ✅ Good: Informative error responses
if gene_id not in self.genes:
    return {'error': f'Gene {gene_id} not found in dataset'}

# ✅ Good: Graceful degradation
if not os.path.exists(self.exp_desc_file):
    return  # Skip loading condition details if file missing

# ❌ Bad: Generic errors
raise Exception("Something went wrong")  # Not helpful for users
```

### Data Return Formats
```python
# ✅ Good: Structured, consistent returns
return {
    'gene': {
        'locusId': gene_info['locusId'],
        'sysName': gene_info['sysName'],
        'description': gene_info['description']
    },
    'fitness_data': fitness_data,
    'total_conditions': len(fitness_data)
}

# ❌ Bad: Inconsistent structure
return gene_info  # Raw internal data structure
```

## File Organization

### Directory Structure
```
fitness-mcp/
├── src/fitness_mcp/
│   └── main.py              # All MCP logic with 3 data loaders
├── data/
│   ├── fit_t.tab            # Main fitness data (90MB TSV)
│   ├── exp_organism_Agro.txt # Experimental condition descriptions
│   ├── RbTnSeq_modules_t1e-7.csv # Gene modules data
│   ├── module_meta.tsv      # Module metadata
│   └── fit_t_pairs_threshold_2_long.tab # Significant gene-condition pairs
├── tests/
│   ├── test_main.py         # Core MCP tool tests
│   └── test_integration.py  # Integration tests
├── prompts/
│   └── fitness-demo-prompt.txt # Claude demo prompt
├── logs/                    # All log files (gitignored)
├── Makefile                 # Build, test, and demo targets
├── pyproject.toml           # Project configuration
├── README.md                # User documentation
└── CLAUDE.md                # Development guide (this file)
```

### Code Organization Within main.py
The monolithic structure contains three main data loaders:

1. **FitnessDataLoader** (lines 214-477) - Main fitness data with thread safety
2. **ModuleDataLoader** (lines 23-212) - Gene functional modules
3. **PairsDataLoader** (lines 563-663) - Significant fitness pairs
4. **MCP Tool Functions** (lines 672-1175) - 16 analysis tools
5. **FastMCP Registration** (lines 1177-1206) - Server setup

All loaders follow TSV best practices with caching, file change detection, and thread safety.

## Testing Standards

### No Mocks or Patches Policy
**NEVER use mocks, patches, or similar test doubles.** Tests should use real implementations and real data.

```python
# ✅ Good: Test with real data loaders and real logic
def test_get_gene_info():
    result = get_gene_info("Atu0001")
    assert result["locusId"] == "Atu0001"

# ❌ Bad: Using mocks or patches
@patch.object(fitness_loader, 'get_gene_info')
def test_get_gene_info_mocked(mock_get):  # DON'T DO THIS
    mock_get.return_value = {'locusId': 'Atu0001'}
```

### Conservative Error Handling in Tests
Be very conservative with try/except in tests. If something fails, we want to know about it immediately.

```python
# ✅ Good: Let exceptions bubble up
def test_gene_lookup():
    result = get_gene_info("Atu0001")  # Will fail clearly if broken
    assert "locusId" in result

# ❌ Bad: Hiding failures with try/except
def test_gene_lookup():
    try:
        result = get_gene_info("Atu0001")
    except Exception:
        assert False  # Doesn't show what actually failed
```

### No Async/Batch Processing
**NEVER use async/await or batch processing.** MCP tools are called by LLM agents (Claude, Goose) that expect synchronous responses.

```python
# ✅ Good: Synchronous tool functions
def get_gene_info(gene_id: str) -> Dict[str, Any]:
    return fitness_loader.get_gene_info(gene_id)

# ❌ Bad: Async functions break MCP protocol
async def get_gene_info_async(gene_id: str) -> Dict[str, Any]:  # DON'T DO THIS
    return await fitness_loader.get_gene_info_async(gene_id)
```

### Test Categories
Use pytest markers for different test types:

```python
@pytest.mark.unit
def test_fitness_score_interpretation():
    """Fast, isolated tests of core logic."""
    
@pytest.mark.integration  
def test_full_gene_analysis():
    """Tests involving multiple components."""
    
@pytest.mark.slow
def test_large_dataset_performance():
    """Performance tests with full 90MB dataset."""
```

### Test Data
```python
# ✅ Good: Small, focused test data
TEST_GENE_DATA = {
    'Atu0001': {
        'locusId': 'Atu0001',
        'sysName': 'test_gene',
        'description': 'test gene description',
        'fitness_values': [-0.5, 0.8, 0.1]
    }
}

# ❌ Bad: Using full production dataset in unit tests
def test_search():
    loader = FitnessDataLoader("data/fit_t.tab")  # 90MB file in tests
```

## Performance Guidelines

### Memory Efficiency
- Load data once, share across all tools
- Use generators for large result sets
- Cache expensive computations with LRU
- Monitor memory usage with large datasets

### Response Time Targets
- Gene lookups: < 10ms
- Simple searches: < 100ms  
- Complex analyses: < 1s
- Data loading: < 5s (one-time cost)

## Security Considerations

### Input Validation
```python
# ✅ Good: Validate all user inputs
if not isinstance(gene_id, str) or not gene_id.strip():
    return {'error': 'Gene ID must be a non-empty string'}

# ✅ Good: Sanitize filter inputs
condition_filter = condition_filter.lower().strip() if condition_filter else None

# ❌ Bad: Direct use of user input in file operations
with open(f"data/{user_input}.tab", 'r'):  # Path injection risk
```

### Data Privacy
- No logging of sensitive biological data
- Error messages should not expose internal data structures
- Rate limiting considerations for MCP deployment

## Dependencies Management

### Core Dependencies
- `fastmcp`: MCP server framework
- `csv`: Built-in TSV parsing (sufficient for current needs)
- `threading`: Thread safety for concurrent access
- `functools.lru_cache`: Performance optimization

### Dependency Principles
- Prefer standard library when possible
- Avoid heavy dependencies (pandas, numpy) unless necessary
- Keep startup time fast
- Consider memory footprint of dependencies

## Documentation Standards

### Docstring Requirements
All public functions must have docstrings following this pattern:

```python
def function_name(param1: str, param2: Optional[float] = None) -> Dict[str, Any]:
    """
    Brief one-line description.

    Longer description explaining biological context and use cases.

    Args:
        param1: Description with biological meaning and examples
        param2: Description with default behavior explanation

    Returns:
        Description of return structure and biological interpretation
        
    Raises:
        SpecificError: When this specific condition occurs
    """
```

### Code Comments
```python
# ✅ Good: Explain biological reasoning
if fitness_val < -0.5:  # Gene knockout improves fitness significantly
    gene_inhibits_growth.append(item)

# ✅ Good: Explain technical decisions  
self._mtime = os.path.getmtime(self.data_file)  # Track file changes for cache invalidation

# ❌ Bad: Obvious comments
gene_id = row[0]  # Get gene ID from first column
```

## Deployment Considerations

### MCP Server Configuration
- Server should start quickly (< 5s including data loading)
- Handle concurrent requests safely
- Graceful error responses for client consumption
- Memory usage should be predictable and bounded

### Production Readiness
- Monitor file system access for data file changes
- Log meaningful events for debugging
- Handle edge cases in biological data gracefully
- Provide clear error messages for common user mistakes

---

## GitHub Issue Workflow

When working on GitHub issues, follow this systematic approach:

1. **Create and checkout linked branch**: `git checkout -b issue-{number}-{short-description}`
2. **Bring feature branch up to date with main**: Feature branches MUST always be brought up to date with main before completing work
3. **Work on implementation** following all guidelines in this document
4. **Test thoroughly** with comprehensive test suite
5. **Run `make all` until all problems are fixed** - This is mandatory before completion
6. **Commit and push** when complete
7. **Create pull request** that auto-closes the issue when merged

### Branch Naming Convention
- Format: `issue-{number}-{short-description}`
- Examples: `issue-21-metadata-registry`, `issue-22-data-redundancy`
- Use kebab-case for descriptions
- Keep descriptions concise but descriptive

### Task Completion Requirements
**ALWAYS run `make all` before considering any task complete.** This command runs the full validation pipeline:
- Code formatting and linting
- Type checking 
- Comprehensive test suite
- Coverage verification
- Integration testing
- MCP protocol validation

**Repeat `make all` until zero errors remain.** Only then is a task ready for commit and push. This ensures:
- No regressions introduced
- All code quality standards met
- Full functionality preserved
- Ready for production deployment

---

## Future Development Guidelines

### Extensibility Patterns
When adding new features:
1. Follow the existing data loader pattern
2. Add appropriate caching for expensive operations
3. Maintain biological accuracy in all interpretations
4. Add comprehensive docstrings with examples
5. Consider thread safety for concurrent access

### Breaking Changes
Any changes to the biological interpretation of fitness scores constitute breaking changes and must be carefully reviewed.

### Code Review Checklist
- [ ] Correct biological interpretation maintained
- [ ] Thread safety considered
- [ ] Appropriate caching implemented
- [ ] Comprehensive docstrings with examples
- [ ] Error handling for edge cases
- [ ] Performance impact assessed
- [ ] Tests added for new functionality