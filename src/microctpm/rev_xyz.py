
# This module was migrated from the original research script: rev_xyz.py
# Hard-coded example execution was removed. Call the functions directly or use examples/.
import numpy as np
import scipy.ndimage as ndi
from skimage import io
import matplotlib.pyplot as plt
import tifffile
import pandas as pd
import os

# Voxel size in µm (corrected to 46×46×46 µm)
VOXEL_SIZE = np.array([15, 15, 15])  # µm
VOXEL_VOLUME = np.prod(VOXEL_SIZE)  # µm³

# 1. Load and filter 3D µCT image
def load_and_filter_image(image_path):
    """Load 3D µCT image and apply Gaussian filter to reduce noise."""
    image = tifffile.imread(image_path)
    if len(image.shape) != 3:
        raise ValueError("Input image must be 3D")
    filtered_image = ndi.gaussian_filter(image, sigma=1.0)
    print(f"Image Intensity Range: min={float(np.min(filtered_image)):.2f}, max={float(np.max(filtered_image)):.2f}")
    print(f"Image Shape: {image.shape}")
    return filtered_image

# 2. Create subvolumes with a fixed threshold of 80 applied to the main image
def create_subvolumes(image, fixed_threshold=85):
    """Create subvolumes from the main image after applying a fixed threshold of 80."""
    if len(image.shape) != 3:
        raise ValueError("Input image must be 3D")
    image_shape = image.shape
    binary_image = image < fixed_threshold  # Apply fixed threshold of 80 to main image (pores < 80 = 1s, solids >= 80 = 0s)
    print(f"Applied fixed threshold of {fixed_threshold} to the full image.")
    print(f"Full Image Porosity after threshold: {np.sum(binary_image) / binary_image.size * 100:.2f}%")
    
    subvolumes = {}
    positions = []
    sub_volume_idx = 0
    division_steps = [2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50]  # Specified division steps
    
    for n in division_steps:
        depth_step = image_shape[0] // n
        height_step = image_shape[1] // n
        width_step = image_shape[2] // n
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    z_start = k * depth_step
                    z_end = (k + 1) * depth_step
                    y_start = i * height_step
                    y_end = (i + 1) * height_step
                    x_start = j * width_step
                    x_end = (j + 1) * width_step
                    z_end = min(z_end, image_shape[0])
                    y_end = min(y_end, image_shape[1])
                    x_end = min(x_end, image_shape[2])
                    sub_binary = binary_image[z_start:z_end, y_start:y_end, x_start:x_end]
                    if sub_binary.shape[0] == 0 or sub_binary.shape[1] == 0 or sub_binary.shape[2] == 0:
                        continue
                    pos_key = f"subvol_n{n}_{sub_volume_idx:03d}_at_z{z_start}_y{y_start}_x{x_start}"
                    subvolumes[pos_key] = {'binary': sub_binary}
                    positions.append(pos_key)
                    porosity = np.sum(sub_binary) / sub_binary.size  # Decimal fraction
                    print(f"Created {pos_key} - Shape: {sub_binary.shape}, Porosity: {porosity:.4f}")
                    sub_volume_idx += 1
    
    return subvolumes, positions

# 3. Compute visible pore volume metrics and plot with CV
def compute_pore_volume_metrics(subvolumes, positions, voxel_size, output_dir):
    """Compute visible pore volume range, mean, std, CV, and plot against subvolume size with CV."""
    size_metrics = {}
    all_subvolume_data = []
    cv_data = []
    unique_sizes = sorted(set(int(pos.split('_')[1][1:]) for pos in positions))
    
    for n in unique_sizes:
        phi_vis_list = []
        subvolume_sizes = []
        for pos in [p for p in positions if f"subvol_n{n}_" in p]:
            sub_binary = subvolumes[pos]['binary']
            if len(sub_binary.shape) != 3:
                raise ValueError(f"Sub-volume {pos} is not 3D")
            dims = sub_binary.shape
            # Use the smallest subdivided dimension as the characteristic size
            size_um = min(dims[0], dims[1], dims[2]) * voxel_size  # Depth, Height, or Width in µm
            phi_vis = np.sum(sub_binary) / sub_binary.size
            phi_vis_list.append(phi_vis)
            subvolume_sizes.append(size_um)
            all_subvolume_data.append({
                'Size (n)': n,
                'Physical Size (µm)': size_um,
                'Phi_vis (fraction)': phi_vis,
                'Position': pos
            })
        
        phi_vis_mean = np.mean(phi_vis_list)
        phi_vis_std = np.std(phi_vis_list)
        phi_vis_cv = phi_vis_std / phi_vis_mean if phi_vis_mean != 0 else float('inf')
        size_metrics[n] = {
            'phi_vis_mean': phi_vis_mean,
            'phi_vis_std': phi_vis_std,
            'phi_vis_cv': phi_vis_cv,
            'phi_vis_values': phi_vis_list,
            'subvolume_sizes': subvolume_sizes
        }
        cv_data.append({
            'Size (n)': n,
            'Mean Physical Size (µm)': np.mean(subvolume_sizes),
            'Mean Phi_vis (fraction)': phi_vis_mean,
            'Std Dev Phi_vis': phi_vis_std,
            'Coefficient of Variation': phi_vis_cv
        })
        print(f"Size n={n}: Mean φ_vis = {phi_vis_mean:.4f}, Std Dev = {phi_vis_std:.4f}, CV = {phi_vis_cv:.3f}")
    
    # Save subvolume data and CV data to CSV files
    pd.DataFrame(all_subvolume_data).to_csv(os.path.join(output_dir, "subvolume_pore_volume_fullxyz_WF4raw_cv_1.csv"), index=False)
    pd.DataFrame(cv_data).to_csv(os.path.join(output_dir, "coefficient_of_variation_WF4raw_1.csv"), index=False)
    
    sizes = list(size_metrics.keys())
    mean_sizes = [np.mean(size_metrics[s]['subvolume_sizes']) for s in sizes]
    phi_vis_means = [size_metrics[s]['phi_vis_mean'] for s in sizes]
    phi_vis_stds = [size_metrics[s]['phi_vis_std'] for s in sizes]
    phi_vis_cvs = [size_metrics[s]['phi_vis_cv'] for s in sizes]
    all_sizes = []
    all_phi_vis = []
    for s in sizes:
        all_sizes.extend(size_metrics[s]['subvolume_sizes'])
        all_phi_vis.extend(size_metrics[s]['phi_vis_values'])
    
    overall_phi_vis_mean = np.mean(all_phi_vis)
    print(f"Overall Mean Visible Pore Volume: {overall_phi_vis_mean:.4f}")
    
    # Plot with improved aesthetics
    fig, ax1 = plt.subplots(figsize=(10, 6))
    
    # Scatter plot for φ_vis
    ax1.scatter(all_sizes, all_phi_vis, color='blue', alpha=0.3, label='Subvolume Data')
    ax1.errorbar(mean_sizes, phi_vis_means, yerr=phi_vis_stds, fmt='b-o', ecolor='black', capsize=5, capthick=2, label='Mean Visible Porosity')
    ax1.axhline(y=overall_phi_vis_mean, color='green', linestyle='--', label=f'Overall Mean φ_vis = {overall_phi_vis_mean:.4f}')
    ax1.set_xlabel('Subvolume Size xyz (µm)')
    ax1.set_ylabel('Visible Pore Volume (fraction)', color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # CV plot
    ax2 = ax1.twinx()
    ax2.plot(mean_sizes, phi_vis_cvs, 'r--o', label='Coefficient of Variation')
    ax2.set_ylabel('Coefficient of Variation', color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.set_ylim(0.0, 0.5)  # Adjust CV axis for better scaling
    
    # Title and legend
    fig.suptitle('Visible Pore Volume and Coefficient of Variation vs. Subvolume Size')
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "pore_volume_and_cv_vs_size_full_PeatS_1_cv_3D.tiff"))
    plt.show()
    
    return size_metrics

# Main execution


def run_rev_xyz(image_path, output_dir, voxel_size=VOXEL_SIZE[0], fixed_threshold=85):
    """Run the full 3D REV workflow from image path to output files."""
    os.makedirs(output_dir, exist_ok=True)
    filtered_image = load_and_filter_image(image_path)
    subvolumes, positions = create_subvolumes(filtered_image, fixed_threshold=fixed_threshold)
    return compute_pore_volume_metrics(subvolumes, positions, voxel_size, output_dir)
