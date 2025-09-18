.PHONY: test-coverage clean install dev format lint all server build upload-test upload release deptry mypy test-fitness-protocol test-gene-analysis test-integration test-version test-gene-fitness test-claude-mcp

# Default target
all: clean install dev test-coverage format lint mypy deptry build test-fitness-protocol test-gene-analysis test-integration test-version

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
	rm -rf logs/

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
test-fitness-protocol:
	@echo "Testing fitness MCP protocol handshake..."
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
	 echo '{"jsonrpc": "2.0", "method": "tools/call", "params": {"name": "get_gene_info", "arguments": {"gene_id": "Atu3150"}}, "id": 2}') | \
	if [ -n "$$TIMEOUT_CMD" ]; then $$TIMEOUT_CMD 5 uv run fitness-mcp; else uv run fitness-mcp & PID=$$!; sleep 5; kill $$PID 2>/dev/null || true; fi

test-gene-analysis:
	@echo "Testing gene fitness analysis via MCP..."
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
	 echo '{"jsonrpc": "2.0", "method": "tools/call", "params": {"name": "search_genes", "arguments": {"query": "ribosome", "limit": 2}}, "id": 3}'; \
	 sleep 0.1; \
	 echo '{"jsonrpc": "2.0", "method": "tools/call", "params": {"name": "find_essential_genes", "arguments": {"condition_filter": "pH", "min_fitness_threshold": 1.0, "limit": 1}}, "id": 4}') | \
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
logs:
	@echo "Creating logs directory..."
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
		2>&1 | tee logs/claude-mcp-test.log

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
		2>&1 | tee logs/atu3150-function-analysis.log
	@echo "✅ Check logs/atu3150-function-analysis.log for biological insights"


# Show all available Claude Code demos
demo-help: ## Show all available Claude Code demo options
	@echo "🧬 FITNESS-MCP CLAUDE CODE DEMO"
	@echo "==============================="
	@echo ""
	@echo "📊 Available demo target:"
	@echo "  make demo-atu3150-function   - Analyze Atu3150 gene function using fitness data"
	@echo ""
	@echo "🎯 Demo creates a log file in logs/ with results for analysis"
	@echo "💡 Shows how to infer gene function from fitness patterns!"


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