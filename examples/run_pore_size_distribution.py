from microctpm.pore_size_distribution import run_pore_size_distribution

run_pore_size_distribution(
    image_path="sample.tif",
    output_dir="outputs/psd_wrc",
    voxel_size=15,
    fixed_threshold=85,
)
