.PHONY: test-coverage clean install dev format lint all server build upload-test upload release deptry mypy test-mcp test-mcp-extended test-integration test-version test-gene-fitness test-biological-analysis demo-claude-code test-claude-mcp

# Default target
all: clean install dev test-coverage format lint mypy deptry build test-mcp test-mcp-extended test-integration test-version

# Install everything for development
dev:
	uv sync --dev

# Install production only
install:
	uv sync

# Run tests with coverage
test-coverage:
	uv run pytest --cov=src/fitness_mcp --cov-report=html --cov-report=term tests/

# Clean up build artifacts
clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf htmlcov/
	rm -f .coverage
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	rm -rf src/*.egg-info

# Run server mode
server:
	uv run fitness-mcp

# Format code with ruff
format:
	uv run ruff format src/ tests/

# Lint code with ruff
lint:
	uv run ruff check --fix src/ tests/

# Check for unused dependencies
deptry:
	uvx deptry .

# Type checking
mypy:
	uv run mypy src/

# Build package with uv
build:
	uv build

# Upload to TestPyPI (using token-based auth - set UV_PUBLISH_TOKEN environment variable first)
upload-test:
	uv publish --publish-url https://test.pypi.org/legacy/

# Upload to PyPI (using token-based auth - set UV_PUBLISH_TOKEN environment variable first)  
upload:
	uv publish

# Complete release workflow
release: clean install test-coverage build

# Integration Testing
test-integration:
	@echo "🧬 Testing fitness MCP integration..."
	uv run pytest tests/test_api.py -v

# Gene Fitness Analysis Testing
test-gene-fitness:
	@echo "🔬 Testing gene fitness analysis..."
	uv run pytest tests/test_main.py::test_gene_fitness_analysis -v -s

# Biological Analysis Testing
test-biological-analysis:
	@echo "🧪 Testing biological insights..."
	uv run pytest tests/test_main.py::test_essential_genes -v -s
	uv run pytest tests/test_main.py::test_growth_inhibitors -v -s

# Demo biological functionality  
demo-biological:
	@echo "🚀 FITNESS-MCP BIOLOGICAL ANALYSIS DEMO"
	@echo "======================================="
	@echo "Finding essential genes in pH conditions..."
	@uv run python -c "from fitness_mcp.main import find_essential_genes; import json; print(json.dumps(find_essential_genes('pH', 0.5, 3), indent=2))" | head -20
	@echo ""
	@echo "Finding growth inhibitor genes in carbon conditions..."
	@uv run python -c "from fitness_mcp.main import find_growth_inhibitor_genes; import json; print(json.dumps(find_growth_inhibitor_genes('carbon', -0.5, 3), indent=2))" | head -20
	@echo ""
	@echo "✅ Fitness-MCP provides structured biological insights for gene function analysis!"

# MCP Server testing
test-mcp:
	@echo "Testing MCP protocol handshake..."
	@if command -v timeout >/dev/null 2>&1; then \
		TIMEOUT_CMD="timeout"; \
	elif command -v gtimeout >/dev/null 2>&1; then \
		TIMEOUT_CMD="gtimeout"; \
	else \
		echo "Warning: timeout command not found, skipping timeout"; \
		TIMEOUT_CMD=""; \
	fi; \
	(echo '{"jsonrpc": "2.0", "method": "initialize", "params": {"protocolVersion": "2025-03-26", "capabilities": {"tools": {}}, "clientInfo": {"name": "test-client", "version": "1.0.0"}}, "id": 1}'; \
	 sleep 0.1; \
	 echo '{"jsonrpc": "2.0", "method": "notifications/initialized", "params": {}}'; \
	 sleep 0.1; \
	 echo '{"jsonrpc": "2.0", "method": "tools/list", "id": 2}') | \
	if [ -n "$$TIMEOUT_CMD" ]; then $$TIMEOUT_CMD 5 uv run fitness-mcp; else uv run fitness-mcp & PID=$$!; sleep 5; kill $$PID 2>/dev/null || true; fi

test-mcp-extended:
	@echo "Testing MCP protocol with fitness analysis..."
	@if command -v timeout >/dev/null 2>&1; then \
		TIMEOUT_CMD="timeout"; \
	elif command -v gtimeout >/dev/null 2>&1; then \
		TIMEOUT_CMD="gtimeout"; \
	else \
		echo "Warning: timeout command not found, skipping timeout"; \
		TIMEOUT_CMD=""; \
	fi; \
	(echo '{"jsonrpc": "2.0", "method": "initialize", "params": {"protocolVersion": "2025-03-26", "capabilities": {"tools": {}}, "clientInfo": {"name": "test-client", "version": "1.0.0"}}, "id": 1}'; \
	 sleep 0.1; \
	 echo '{"jsonrpc": "2.0", "method": "notifications/initialized", "params": {}}'; \
	 sleep 0.1; \
	 echo '{"jsonrpc": "2.0", "method": "tools/call", "params": {"name": "search_genes", "arguments": {"query": "ribosome", "limit": 3}}, "id": 3}'; \
	 sleep 0.1; \
	 echo '{"jsonrpc": "2.0", "method": "tools/call", "params": {"name": "find_essential_genes", "arguments": {"condition_filter": "pH", "min_fitness_threshold": 0.5, "limit": 2}}, "id": 4}') | \
	if [ -n "$$TIMEOUT_CMD" ]; then $$TIMEOUT_CMD 10 uv run fitness-mcp; else uv run fitness-mcp & PID=$$!; sleep 10; kill $$PID 2>/dev/null || true; fi

# Test version flag
test-version:
	@echo "🔢 Testing package installation..."
	@echo "Package version:"
	@uv run python -c "import fitness_mcp; print('fitness_mcp package loaded successfully')"

# Data validation tests
test-data:
	@echo "📊 Testing fitness data loading..."
	@uv run python -c "from fitness_mcp.main import fitness_loader; fitness_loader.load_data(); print(f'Loaded {len(fitness_loader.genes)} genes across {len(fitness_loader.conditions)} conditions')"

# Performance testing with full dataset
test-performance:
	@echo "⚡ Testing performance with full dataset..."
	@time uv run python -c "from fitness_mcp.main import fitness_loader, search_genes; fitness_loader.load_data(); print('Data loaded'); results = search_genes('ribosome', 10); print(f'Found {len(results)} genes')"

# Directory creation targets (not .PHONY since they create directories)
data/outputs:
	@echo "Creating outputs directory..."
	@mkdir -p $@

prompts:
	@echo "Creating prompts directory..."
	@mkdir -p $@

# Test with Claude CLI using local MCP config
test-claude-mcp:
	@echo "🧬 Testing fitness-mcp with Claude CLI..."
	claude \
		--debug \
		--verbose \
		--mcp-config .mcp.json \
		--dangerously-skip-permissions \
		--print "Test the fitness-mcp by listing available tools and then search for genes containing 'ribosome' and return the top 3 results" \
		2>&1 | tee claude-mcp-test.log

# Analyze specific gene function using fitness data (Atu3150 lactose transporter example)
demo-atu3150-function: prompts/fitness-demo-prompt.txt ## Analyze what Atu3150 does using fitness data
	@echo "🧬 Analyzing Atu3150 gene function using fitness data..."
	@echo "This demo shows how fitness analysis reveals Atu3150 is a lactose transporter"
	claude \
		--debug \
		--verbose \
		--mcp-config .mcp.json \
		--dangerously-skip-permissions \
		--print "$(shell cat prompts/fitness-demo-prompt.txt)" \
		2>&1 | tee atu3150-function-analysis.log
	@echo "✅ Check atu3150-function-analysis.log for biological insights"

# Comprehensive biological analysis showcasing all fitness MCP capabilities
demo-biological-analysis: prompts/comprehensive-fitness-demo.txt ## Comprehensive biological gene analysis
	@echo "🧬 Running comprehensive biological analysis demo..."
	@echo "This demo showcases:"
	@echo "  - Gene function analysis through fitness patterns"
	@echo "  - Essential gene identification in specific conditions"
	@echo "  - Gene-condition network relationships" 
	@echo "  - Functional module analysis"
	@echo "  - Cross-referencing annotations with experimental data"
	claude \
		--debug \
		--verbose \
		--mcp-config .mcp.json \
		--dangerously-skip-permissions \
		--output-format json \
		--print "$(shell cat prompts/comprehensive-fitness-demo.txt)" \
		2>&1 | tee biological-analysis-demo.log
	@echo "✅ Check biological-analysis-demo.log for comprehensive analysis results"

# Gene-condition network expansion and interaction analysis
demo-network-expansion: prompts/network-analysis-demo.txt ## Analyze gene-condition interaction networks
	@echo "🔗 Running gene-condition network expansion demo..."
	@echo "This demo showcases:"
	@echo "  - Gene-condition fitness relationships"
	@echo "  - Two-hop network expansion (genes→conditions→genes)"
	@echo "  - Biological network interpretation"
	@echo "  - Functional relationship discovery"
	claude \
		--debug \
		--verbose \
		--mcp-config .mcp.json \
		--dangerously-skip-permissions \
		--output-format json \
		--print "$(shell cat prompts/network-analysis-demo.txt)" \
		2>&1 | tee network-expansion-demo.log
	@echo "✅ Check network-expansion-demo.log for network analysis results"

# Show all available Claude Code demos
demo-help: ## Show all available Claude Code demo options
	@echo "🧬 FITNESS-MCP CLAUDE CODE DEMOS"
	@echo "================================"
	@echo ""
	@echo "📊 Available demo targets:"
	@echo "  make demo-atu3150-function   - Analyze Atu3150 gene function using fitness data"
	@echo "  make demo-biological-analysis - Comprehensive biological gene analysis"
	@echo "  make demo-network-expansion  - Analyze gene-condition interaction networks"
	@echo "  make demo-claude-mcp         - Test MCP protocol integration"
	@echo ""
	@echo "🎯 Each demo creates a log file with results for analysis"
	@echo "💡 Use --output-format json for structured data output"

# Backward compatibility aliases (deprecated - use new names above)
demo-claude-code: demo-atu3150-function
demo-claude-comprehensive: demo-biological-analysis  
demo-claude-network: demo-network-expansion

# FITNESS MCP - Claude Desktop config:
#   Add to ~/Library/Application Support/Claude/claude_desktop_config.json:
#   {
#     "mcpServers": {
#       "fitness-mcp": {
#         "command": "uvx",
#         "args": ["fitness-mcp"]
#       }
#     }
#   }
#
# Claude Code MCP setup:
#   claude mcp add -s project fitness-mcp uvx fitness-mcp
#
# Goose setup:
#   goose session --with-extension "uvx fitness-mcp"