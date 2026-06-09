from microctpm.rev_xyz import run_rev_xyz

results = run_rev_xyz(
    image_path="sample.tif",
    output_dir="outputs/rev_xyz",
    voxel_size=15,
    fixed_threshold=85,
)

print("REV XYZ completed.")
print(results.keys())
