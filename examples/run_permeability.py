from microctpm.permeability import calculate_permeability

image_path = "sample.tif"
output_dir = "outputs/permeability"

result = calculate_permeability(
    image_path=image_path,
    output_dir=output_dir,
    voxel_size=15e-6,
    threshold=85,
    inlet_pressure=2000.0,
    outlet_pressure=1990.0,

    # For this validation image, 1990 Pa reproduces the original research-script
    # pressure-drop convention. Users should replace this with their own ΔP.
    pressure_drop=1990.0,

    save_outputs=True,
    make_plot=True,
    show_plot=False,
)

print(result)
