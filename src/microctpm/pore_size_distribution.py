# Migrated from original pore_size_distribution.py
import numpy as np
import scipy.ndimage as ndi
import porespy as ps
import tifffile
import pandas as pd
import os

VOXEL_SIZE = np.array([15, 15, 15], dtype=np.float64)  # Voxel size in µm
VOXEL_VOLUME = np.prod(VOXEL_SIZE)  # = 46 * 46 * 46 = 97,336 µm³
# Absolute path for saving files
OUTPUT_DIR = "results"
os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_and_filter_image(image_path):
    """Load a 3D µCT image and apply Gaussian filter."""
    image = tifffile.imread(image_path)
    filtered_image = ndi.gaussian_filter(image, sigma=1)
    print(
        f"Image at {image_path} loaded and filtered. Shape: {filtered_image.shape}")
    return filtered_image


def compute_pore_volumes_and_sizes(image, output_dir=None, voxel_size=None, fixed_threshold=85, output_filename="pore_size_distribution_wrc.xlsx"):
    global VOXEL_SIZE, VOXEL_VOLUME, OUTPUT_DIR
    if voxel_size is not None:
        VOXEL_SIZE = np.array([voxel_size, voxel_size, voxel_size], dtype=np.float64) if np.isscalar(voxel_size) else np.array(voxel_size, dtype=np.float64)
        VOXEL_VOLUME = np.prod(VOXEL_SIZE)
    if output_dir is not None:
        OUTPUT_DIR = str(output_dir)
        os.makedirs(OUTPUT_DIR, exist_ok=True)
    """Compute total pore volume, pore size distribution, and water retention curve using PoreSpy."""
    # Apply fixed threshold of 75
    print(f"Using Fixed Threshold: {fixed_threshold}")
    binary_image = image < fixed_threshold  # True for pores, False for solids

    # Calculate total pore volume
    pore_voxel_count = np.sum(binary_image)
    total_pore_volume = pore_voxel_count * VOXEL_VOLUME  # in µm³
    total_sample_volume = binary_image.size * VOXEL_VOLUME  # in µm³

    # Compute pore size distribution using PoreSpy
    lt = ps.filters.local_thickness(binary_image)
    psd = ps.metrics.pore_size_distribution(
        lt, bins=100, log=False, voxel_size=np.mean(VOXEL_SIZE))

    # Extract pore sizes (diameters) and volumes
    single_pore_radii = psd.bin_centers  # Radii in µm (scaled by voxel_size)
    single_pore_diameters = single_pore_radii * 2  # Convert to diameters
    # Calculate pore volumes from pdf and total pore volume
    single_pore_volumes = psd.pdf * total_pore_volume * \
        psd.bin_widths  # Volume in each bin (µm³)

    # Filter out zero or invalid volumes
    valid_mask = single_pore_volumes > 0
    single_pore_diameters = single_pore_diameters[valid_mask]
    single_pore_volumes = single_pore_volumes[valid_mask]

    # Prepare metrics
    metrics = {
        'Total Pore Volume (µm³)': float(total_pore_volume),
        'Total Sample Volume (µm³)': float(total_sample_volume),
        'Number of Pores': len(single_pore_volumes),
        'Average Single Pore Volume (µm³)': float(np.mean(single_pore_volumes)) if single_pore_volumes.size > 0 else 0.0,
        'Average Single Pore Size (µm)': float(np.mean(single_pore_diameters)) if single_pore_diameters.size > 0 else 0.0
    }

    # Save single pore volumes and sizes to DataFrame
    pore_data = pd.DataFrame({
        'Single Pore Volume (µm³)': single_pore_volumes,
        'Single Pore Size (µm)': single_pore_diameters
    })

    # Compute water retention curve using simplified Young-Laplace equation
    if single_pore_diameters.size > 0:
        # Calculate matric potential using psi = 2980 / d (psi in cm, d in µm)
        matric_potentials_cm = 2980 / single_pore_diameters  # in cm of water

        # Sort pores by size (ascending) and corresponding volumes
        sorted_indices = np.argsort(single_pore_diameters)
        sorted_pore_diameters = single_pore_diameters[sorted_indices]
        sorted_pore_volumes = single_pore_volumes[sorted_indices]
        sorted_matric_potentials_cm = matric_potentials_cm[sorted_indices]

        # Compute cumulative pore volume for water retention curve
        cumulative_volumes = np.cumsum(
            sorted_pore_volumes)  # Cumulative volume in µm³
        volumetric_water_content = cumulative_volumes / total_sample_volume  # θ

        # Create WRC DataFrame
        wrc_data = pd.DataFrame({
            'Pore Size (µm)': sorted_pore_diameters,
            'Matric Potential (cm H2O)': sorted_matric_potentials_cm,
            'Volumetric Water Content': volumetric_water_content
        })
    else:
        print("No pores detected, skipping WRC calculation.")
        wrc_data = pd.DataFrame({
            'Pore Size (µm)': [],
            'Matric Potential (cm H2O)': [],
            'Volumetric Water Content': []
        })

    # Save to Excel
    if OUTPUT_DIR is None:
        OUTPUT_DIR = os.getcwd()
    output_path = os.path.join(OUTPUT_DIR, output_filename)
    with pd.ExcelWriter(output_path) as writer:
        pd.DataFrame([metrics]).to_excel(
            writer, sheet_name="Summary", index=False)
        pore_data.to_excel(writer, sheet_name="Pore Data", index=False)
        wrc_data.to_excel(
            writer, sheet_name="Water Retention Curve", index=False)
    print(f"Saved {output_path}")

    print("\nPore Metrics:")
    for key, value in metrics.items():
        print(f"{key}: {float(value):.6f}" if isinstance(
            value, float) else f"{key}: {value}")

    return metrics, pore_data, wrc_data


def run_pore_size_distribution(image_path, output_dir, voxel_size=46, fixed_threshold=85, output_filename="pore_size_distribution_wrc.xlsx"):
    """Run the full pore-size-distribution and water-retention workflow."""
    filtered_image = load_and_filter_image(image_path)
    return compute_pore_volumes_and_sizes(filtered_image, output_dir=output_dir, voxel_size=voxel_size, fixed_threshold=fixed_threshold, output_filename=output_filename)
