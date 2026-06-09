from microctpm.critical_pore import calculate_critical_pore_diameter

critical_diameter = calculate_critical_pore_diameter(
    image_path="sample.tif",
    voxel_size=15e-6,
    threshold=85,
    r_max=10,
    sigma=1,
    output_dir="outputs/critical_pore",
    save_plot=True,
    show_plot=False,
)

print("Critical pore diameter:", critical_diameter, "µm")
