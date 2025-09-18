"""Shared test fixtures for fitness_mcp tests."""

import pytest
import os
from src.fitness_mcp.main import metadata_registry


@pytest.fixture(scope="session")
def loaded_metadata_registry():
    """Load metadata registry once per test session to speed up tests.

    This fixture loads the real data files once and reuses the loaded registry
    across all tests in the session. If data files are not available,
    tests that depend on this fixture will be skipped.
    """
    # Check if data files exist
    data_files_exist = all(
        [
            os.path.exists("data/fit_t.tab"),
            os.path.exists("data/exp_organism_Agro.txt"),
            os.path.exists("data/RbTnSeq_modules_t1e-7.csv"),
            os.path.exists("data/module_meta.tsv"),
        ]
    )

    if not data_files_exist:
        pytest.skip("Real fitness data files not available for session")

    # Load data once for the entire session
    metadata_registry.load_data()

    # Verify we have loaded data
    if not metadata_registry.loaded or len(metadata_registry.genes) == 0:
        pytest.skip("Failed to load fitness data for session")

    return metadata_registry


@pytest.fixture(scope="session")
def sample_gene_id(loaded_metadata_registry):
    """Get a sample gene ID from the loaded data for tests."""
    return next(iter(loaded_metadata_registry.genes.keys()))


@pytest.fixture(scope="session")
def sample_condition_id(loaded_metadata_registry):
    """Get a sample condition ID from the loaded data for tests."""
    return next(iter(loaded_metadata_registry.conditions.keys()))


@pytest.fixture(scope="session")
def sample_module_id(loaded_metadata_registry):
    """Get a sample module ID from the loaded data for tests."""
    if loaded_metadata_registry.modules:
        return next(iter(loaded_metadata_registry.modules.keys()))
    return None
