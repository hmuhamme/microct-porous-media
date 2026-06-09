
def test_import_microctpm():
    import microctpm
    assert microctpm.__version__ == "1.0.0"

def test_import_functions():
    from microctpm.rev_xy import run_rev_xy
    from microctpm.rev_xyz import run_rev_xyz
    from microctpm.critical_pore import calculate_critical_pore_diameter
    from microctpm.pore_size_distribution import compute_pore_volumes_and_sizes
    from microctpm.permeability import calculate_permeability
