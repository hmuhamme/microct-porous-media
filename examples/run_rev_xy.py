from microctpm.rev_xy import run_rev_xy

results = run_rev_xy(
    image_path="sample.tif",
    output_dir="outputs/rev_xy",
    voxel_size=15,
    fixed_threshold=85,
)

print("REV XY completed.")
print(results.keys())
