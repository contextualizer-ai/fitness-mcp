# fitness_mcp

A fitness MCP

## Installation

Install the package with uv:

```bash
uv sync --dev
```

## Usage

You can run the MCP:

```bash
uv run fitness-mcp
```

Or install globally and run:

```bash
uvx fitness-mcp
```

Or import in your Python code:

```python
from fitness_mcp.main import create_mcp

mcp = create_mcp()
mcp.run()
```

## Development

### Local Setup

```bash
# Clone the repository
git clone https://github.com/justaddcoffee/fitness-mcp.git
cd fitness-mcp

# Install development dependencies
uv sync --dev
```


## License

BSD-3-Clause
