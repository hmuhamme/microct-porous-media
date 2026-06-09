# Migrated from original pore_metrics.py
import numpy as np
import scipy.ndimage as ndi
from skimage import io, measure, morphology
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import tifffile
import pandas as pd
import os
import porespy as ps
import openpnm as op
from skan import Skeleton, summarize

VOXEL_SIZE = np.array([15, 15, 15], dtype=np.float64)  # Voxel size in µm
VOXEL_VOLUME = np.prod(VOXEL_SIZE)
FIXED_THRESHOLD = 85  # Configurable through compute_metrics/run_pore_metrics
OUTPUT_DIR = "results"  # Default relative output folder
os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_and_filter_image(image_path):
    """Load a single 3D µCT image and apply Gaussian filter to reduce noise."""
    try:
        image = tifffile.imread(image_path)
        if len(image.shape) != 3:
            raise ValueError(f"Input image at {image_path} must be 3D")
        filtered_image = ndi.gaussian_filter(image, sigma=2)
        print(f"Image at {image_path} Intensity Range: min={float(np.min(filtered_image)):.2f}, max={float(np.max(filtered_image)):.2f}")
        print(f"Image Shape: {image.shape}")
        return filtered_image
    except Exception as e:
        print(f"Error loading or filtering image: {e}")
        raise

def compute_euler_number(binary_image):
    try:
        labeled_clusters = measure.label(binary_image, connectivity=3)
        N = np.max(labeled_clusters)
        if N == 0:
            N = len(np.unique(labeled_clusters)) - 1
        skeleton = morphology.skeletonize(binary_image)
        C = np.sum(skeleton > 0) - N
        chi = N - C
        volume = binary_image.size * VOXEL_VOLUME / 1e9
        chi_density = chi / volume if volume > 0 else 0
        print(f"Debug: chi = {chi}, volume = {volume} mm³, chi_density = {chi_density} mm⁻³")
        return chi_density
    except Exception as e:
        print(f"Error computing Euler number: {e}")
        return 0

def compute_gamma_indicator(binary_image):
    try:
        labeled_clusters = measure.label(binary_image)
        N_p = np.sum(binary_image, dtype=np.int64)
        if N_p == 0:
            print("Warning: No pore voxels detected, gamma set to 0")
            return 0
        cluster_sizes = np.bincount(labeled_clusters.ravel())[1:]
        gamma = np.sum(cluster_sizes**2, dtype=np.float64) / (N_p**2)
        return gamma
    except Exception as e:
        print(f"Error computing gamma indicator: {e}")
        return 0

def analyze_pore_network(binary_image):
    """Extract pore network using PoreSpy and OpenPNM."""
    if len(binary_image.shape) != 3:
        raise ValueError("Binary image must be 3D")
    
    binary_image = binary_image.astype(bool)
    
    try:
        print("Extracting pore network with PoreSpy snow2 algorithm...")
        net = ps.networks.snow2(binary_image, voxel_size=VOXEL_SIZE[0] / 1e6)
        print("Snow2 network keys:", list(net.network.keys()))
        pore_diameter = net.network.get('pore.diameter', 
                        net.network.get('pore.equivalent_diameter', 
                        np.zeros(net.network['pore.coords'].shape[0])))
        throat_diameter = net.network.get('throat.diameter', 
                          net.network.get('throat.equivalent_diameter', 
                          np.zeros(net.network['throat.conns'].shape[0])))
        if np.max(pore_diameter) > 1000:  # Assume µm if > 1000
            pore_diameter /= 1e6
            throat_diameter /= 1e6
        net_dict = {
            'pore.coords': net.network['pore.coords'],
            'throat.conns': net.network['throat.conns'],
            'pore.equivalent_diameter': pore_diameter,
            'throat.equivalent_diameter': throat_diameter,
        }
    except Exception as e:
        print(f"Error in PoreSpy snow2 extraction: {e}. Falling back to basic metrics.")
        skeleton = morphology.skeletonize(binary_image)
        num_branches = np.sum(skeleton > 0)
        print(f"Fallback: Number of Skeleton Voxels (Branches): {num_branches}")
        return (num_branches, 0, 0, 1.0, np.zeros(50), np.linspace(0, 2500, 50), 0, [], [])

    try:
        proj = op.Project()
        try:
            pn = op.network.Network(name='network', project=proj)
        except AttributeError:
            try:
                pn = op.network.NetworkBase(name='network', project=proj)
            except AttributeError:
                pn = op.network.GenericNetwork(name='network', project=proj)
        pn.update(net_dict)
        print("OpenPNM network created successfully.")
        print(pn)
    except Exception as e:
        print(f"Error in OpenPNM network creation: {e}. Falling back to basic metrics.")
        skeleton = morphology.skeletonize(binary_image)
        num_branches = np.sum(skeleton > 0)
        print(f"Fallback: Number of Skeleton Voxels (Branches): {num_branches}")
        return (num_branches, 0, 0, 1.0, np.zeros(50), np.linspace(0, 2500, 50), 0, [], [])

    num_pores = pn.Np
    num_throats = pn.Nt
    avg_pore_diameter = np.mean(pn['pore.equivalent_diameter']) * 1e6 if 'pore.equivalent_diameter' in pn and np.any(pn['pore.equivalent_diameter']) else 0
    avg_throat_diameter = np.mean(pn['throat.equivalent_diameter']) * 1e6 if 'throat.equivalent_diameter' in pn and np.any(pn['throat.equivalent_diameter']) else 0
    coordination_number = np.mean(pn.num_neighbors(pn.pores())) if num_pores > 0 else 0

    throat_radii = pn['throat.equivalent_diameter'] * 1e6 / 2 if 'throat.equivalent_diameter' in pn and np.any(pn['throat.equivalent_diameter']) else np.array([])
    if len(throat_radii) > 0:
        hist, bins = np.histogram(throat_radii, bins=50, range=(0, 2500))
        throat_size_freq = hist
        throat_bin_centers = (bins[:-1] + bins[1:]) / 2
    else:
        throat_size_freq = np.zeros(50)
        throat_bin_centers = np.linspace(0, 2500, 50)
        print("Warning: No throat radii detected.")

    plt.bar(throat_bin_centers, throat_size_freq, width=np.diff(throat_bin_centers)[0], align='center')
    plt.title("Pore-Throat Size Distribution")
    plt.xlabel("Throat Diameter (µm)")
    plt.ylabel("Frequency")
    output_path = os.path.join(OUTPUT_DIR, "throat_size_distribution.png")
    plt.savefig(output_path)
    plt.close()
    if not os.path.exists(output_path):
        print(f"Error: Failed to save {output_path}")

    throat_data = pd.DataFrame({'Throat Diameter (µm)': throat_radii})
    output_path = os.path.join(OUTPUT_DIR, "throat_radii.csv")
    try:
        throat_data.to_csv(output_path, index=False)
        if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            print(f"Saved {output_path}")
        else:
            print(f"Error: Failed to save {output_path} or file is empty")
    except Exception as e:
        print(f"Error saving {output_path}: {e}")

    tortuosity_values = []
    if tortuosity_values:
        plt.hist(tortuosity_values, bins=20, color='purple', alpha=0.7)
        plt.title("Tortuosity Distribution")
        plt.xlabel("Tortuosity")
        plt.ylabel("Frequency")
        output_path = os.path.join(OUTPUT_DIR, "tortuosity_distribution.png")
        plt.savefig(output_path)
        plt.close()
        if not os.path.exists(output_path):
            print(f"Error: Failed to save {output_path}")
        tortuosity_data = pd.DataFrame({'Tortuosity': tortuosity_values})
        output_path = os.path.join(OUTPUT_DIR, "tortuosity_values.csv")
        try:
            tortuosity_data.to_csv(output_path, index=False)
            if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
                print(f"Saved {output_path}")
            else:
                print(f"Error: Failed to save {output_path} or file is empty")
        except Exception as e:
            print(f"Error saving {output_path}: {e}")

    print(f"Pore Network Characteristics:")
    print(f"Number of Pores: {num_pores}")
    print(f"Number of Throats: {num_throats}")
    print(f"Average Pore Diameter: {float(avg_pore_diameter):.2f} µm")
    print(f"Average Throat Diameter: {float(avg_throat_diameter):.2f} µm")
    print(f"Average Coordination Number: {float(coordination_number):.2f}")

    return (num_pores, num_throats, avg_pore_diameter, coordination_number,
            throat_size_freq, throat_bin_centers, avg_throat_diameter, throat_radii, tortuosity_values)

def analyze_pore_properties(binary_image):
    if len(binary_image.shape) != 3:
        raise ValueError("Binary image must be 3D")
    print(f"Analyzing full image of shape {binary_image.shape} for pore analysis.")
    total_voxels = binary_image.size
    sub_volume = total_voxels * VOXEL_VOLUME / 1e6  # Convert to mm³
    pore_voxel_count = np.sum(binary_image, dtype=np.int64)
    pore_volume = pore_voxel_count * VOXEL_VOLUME
    pore_percentage = pore_voxel_count / total_voxels * 100 if total_voxels > 0 else 0
    if pore_voxel_count == 0:
        print("Warning: No pore voxels detected.")
        return (0, 0, 0, 0, 0, 0, 0, False, 0, np.zeros(50), np.linspace(0, 2500, 50), 0, [], [], [], 0, 0, 0)
    labeled_pores = measure.label(binary_image)
    regions = measure.regionprops(labeled_pores)
    total_pore_area = 0
    max_diameter = 0
    min_diameter = float('inf')
    pore_shape_factors = []
    pore_diameters = []
    for region in regions:
        boundary = morphology.binary_erosion(binary_image) ^ binary_image
        boundary_region = boundary * (labeled_pores == region.label)
        surface_voxels = np.sum(boundary_region, dtype=np.int64)
        voxel_surface_area = 6 * (VOXEL_SIZE[0] ** 2) / 1e6  # Convert µm² to mm²
        region_surface_area = surface_voxels * voxel_surface_area
        total_pore_area += region_surface_area
        equiv_diameter = region.equivalent_diameter * np.mean(VOXEL_SIZE)
        pore_diameters.append(equiv_diameter)
        max_diameter = max(max_diameter, equiv_diameter)
        min_diameter = min(min_diameter, equiv_diameter) if equiv_diameter > 0 else min_diameter
        volume = region.area * VOXEL_VOLUME / 1e9  # Convert µm³ to mm³
        surface_area = region_surface_area
        if volume > 0 and surface_area > 0:
            sphericity = (np.pi ** (1/3) * (6 * volume) ** (2/3)) / surface_area
            sphericity = min(sphericity, 1.0)
            pore_shape_factors.append(sphericity)
        else:
            pore_shape_factors.append(0)
    min_diameter = min_diameter if min_diameter != float('inf') else 0
    avg_pore_sphericity = np.mean(pore_shape_factors) if pore_shape_factors else 0
    pore_surface_density = total_pore_area / sub_volume if sub_volume > 0 else 0
    distance_map = ndi.distance_transform_edt(binary_image)
    pore_sizes = distance_map[binary_image] * 2 * np.mean(VOXEL_SIZE)
    if len(pore_sizes) == 0:
        print("Warning: Pore sizes empty despite pore voxels. Returning defaults.")
        return (pore_volume, total_pore_area, max_diameter, min_diameter, avg_pore_sphericity,
                pore_surface_density, 0, False, 0, np.zeros(50), np.linspace(0, 2500, 50), 0, [], pore_diameters, pore_shape_factors, 0, 0, 0)
    gradient = np.gradient(distance_map)
    norm_gradient = np.sqrt(sum(g**2 for g in gradient))
    norm_gradient[norm_gradient == 0] = 1e-6
    unit_normals = [g / norm_gradient for g in gradient]
    divergence = sum(np.gradient(un, axis=i) for i, un in enumerate(unit_normals))
    mean_curvature = -0.5 * divergence
    mean_curvature_density = np.mean(mean_curvature[binary_image]) * VOXEL_SIZE[0] / 1e6 if sub_volume > 0 else 0
    labeled, num_clusters = measure.label(binary_image, return_num=True)
    percolates = False
    for cluster_id in range(1, num_clusters + 1):
        cluster = (labeled == cluster_id)
        if np.any(cluster[0, :, :]) and np.any(cluster[-1, :, :]):
            percolates = True
            break
        if np.any(cluster[:, 0, :]) and np.any(cluster[:, -1, :]):
            percolates = True
            break
        if np.any(cluster[:, :, 0]) and np.any(cluster[:, :, -1]):
            percolates = True
            break
    sorted_diameters = np.sort(np.unique(pore_sizes))
    critical_diameter = None
    for d in sorted_diameters:
        filtered_pores = distance_map * (distance_map >= d / (2 * np.mean(VOXEL_SIZE)))
        filtered_binary = filtered_pores > 0
        labeled_filt, num_filt = measure.label(filtered_binary, return_num=True)
        percolates_filt = False
        for cluster_id in range(1, num_filt + 1):
            cluster = (labeled_filt == cluster_id)
            if (np.any(cluster[0, :, :]) and np.any(cluster[-1, :, :])) or \
               (np.any(cluster[:, 0, :]) and np.any(cluster[:, -1, :])) or \
               (np.any(cluster[:, :, 0]) and np.any(cluster[:, :, -1])):
                percolates_filt = True
                break
        if percolates_filt:
            critical_diameter = d
            break
    critical_diameter = critical_diameter if critical_diameter is not None else max_diameter
    hist, bins = np.histogram(pore_sizes, bins=50, range=(0, 2500))
    bin_centers = (bins[:-1] + bins[1:]) / 2
    pore_size_freq = hist
    avg_diameter = np.mean(pore_sizes) if len(pore_sizes) > 0 else 0
    if pore_diameters:
        plt.hist(pore_diameters, bins=20, color='blue', alpha=0.7)
        plt.title("Pore Diameter Distribution")
        plt.xlabel("Pore Diameter (µm)")
        plt.ylabel("Frequency")
        output_path = os.path.join(OUTPUT_DIR, "pore_diameter_distribution.png")
        plt.savefig(output_path)
        plt.close()
        if not os.path.exists(output_path):
            print(f"Error: Failed to save {output_path}")
    if pore_shape_factors:
        plt.hist(pore_shape_factors, bins=20, color='green', alpha=0.7)
        plt.title("Pore Sphericity Distribution")
        plt.xlabel("Sphericity")
        plt.ylabel("Frequency")
        output_path = os.path.join(OUTPUT_DIR, "pore_sphericity_distribution.png")
        plt.savefig(output_path)
        plt.close()
        if not os.path.exists(output_path):
            print(f"Error: Failed to save {output_path}")
    print(f"Debug: Number of pore diameters: {len(pore_diameters)}")
    print(f"Debug: Number of pore sizes: {len(pore_sizes)}")
    pore_data_raw = pd.DataFrame({
        'Pore Diameter (µm)': pore_diameters,
        'Sphericity': pore_shape_factors
    })
    output_path = os.path.join(OUTPUT_DIR, "pore_properties_raw.csv")
    try:
        pore_data_raw.to_csv(output_path, index=False)
        if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            print(f"Saved {output_path}")
        else:
            print(f"Error: Failed to save {output_path} or file is empty")
    except Exception as e:
        print(f"Error saving {output_path}: {e}")
    if len(pore_sizes) == 0:
        print("Error: pore_sizes array is empty, cannot save pore_sizes_raw.csv")
    else:
        # Save pore_sizes in chunks to handle large arrays
        chunk_size = 1000000
        output_path = os.path.join(OUTPUT_DIR, "pore_sizes_raw.csv")
        try:
            with open(output_path, 'w') as f:
                f.write('Pore Size (µm)\n')
                for i in range(0, len(pore_sizes), chunk_size):
                    chunk = pore_sizes[i:i + chunk_size]
                    np.savetxt(f, chunk, fmt='%.6f')
            if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
                print(f"Saved {output_path}")
            else:
                print(f"Error: Failed to save {output_path} or file is empty")
        except Exception as e:
            print(f"Error saving {output_path}: {e}")
    # Compute macro, meso, micro volumes
    macro_pores = pore_sizes > 300
    meso_pores = (pore_sizes >= 50) & (pore_sizes <= 300)
    micro_pores = pore_sizes < 50
    num_macro = np.sum(macro_pores, dtype=np.int64)
    num_meso = np.sum(meso_pores, dtype=np.int64)
    num_micro = np.sum(micro_pores, dtype=np.int64)
    macro_volume = num_macro * VOXEL_VOLUME
    meso_volume = num_meso * VOXEL_VOLUME
    micro_volume = num_micro * VOXEL_VOLUME
    print(f"Pore Properties:")
    print(f"Total Pore Volume: {float(pore_volume):.2f} µm³")
    print(f"Total Pore Area: {float(total_pore_area):.2f} mm²")
    print(f"Maximum Pore Diameter: {float(max_diameter):.2f} µm")
    print(f"Minimum Pore Diameter: {float(min_diameter):.2f} µm")
    print(f"Average Sphericity: {float(avg_pore_sphericity):.2f}")
    print(f"Pore Surface Density: {float(pore_surface_density):.6f} mm⁻¹")
    print(f"Mean Curvature Density: {float(mean_curvature_density):.6f} mm⁻²")
    print(f"Percolation: {'Yes' if percolates else 'No'}")
    print(f"Critical Pore Diameter: {float(critical_diameter):.2f} µm")
    print(f"Average Pore Diameter: {float(avg_diameter):.2f} µm")
    return (pore_volume, total_pore_area, max_diameter, min_diameter, avg_pore_sphericity,
            pore_surface_density, mean_curvature_density, percolates, critical_diameter,
            pore_size_freq, bin_centers, avg_diameter, pore_sizes, pore_diameters, pore_shape_factors,
            macro_volume, meso_volume, micro_volume)

def analyze_pore_structure(binary_image, original_image):
    if len(binary_image.shape) != 3 or len(original_image.shape) != 3:
        raise ValueError("Binary or original image must be 3D")
    print(f"Analyzing full image of shape {binary_image.shape} for pore structure.")
    print(f"Fixed Threshold: {FIXED_THRESHOLD}")
    (pore_volume, total_pore_area, max_diameter, min_diameter, avg_pore_sphericity,
     pore_surface_density, mean_curvature_density, percolates, critical_diameter,
     pore_size_freq, bin_centers, avg_diameter, pore_sizes, pore_diameters, pore_shape_factors,
     macro_volume, meso_volume, micro_volume) = analyze_pore_properties(binary_image)
    print(f"Debug: Total pore voxels: {np.sum(binary_image, dtype=np.int64)}")
    print(f"Debug: Macropore voxels (>300 µm): {np.sum(pore_sizes > 300, dtype=np.int64)}")
    print(f"Debug: Mesopore voxels (50–300 µm): {np.sum((pore_sizes >= 50) & (pore_sizes <= 300), dtype=np.int64)}")
    print(f"Debug: Micropore voxels (<50 µm): {np.sum(pore_sizes < 50, dtype=np.int64)}")
    print(f"Debug: Sum of categorized voxels: {np.sum(pore_sizes > 0, dtype=np.int64)}")
    total_pore_volume = pore_volume
    if total_pore_volume > 0:
        macro_pct = (macro_volume / total_pore_volume) * 100
        meso_pct = (meso_volume / total_pore_volume) * 100
        micro_pct = (micro_volume / total_pore_volume) * 100
    else:
        macro_pct = meso_pct = micro_pct = 0.0
    print(f"Macropores (>300 µm): {float(macro_volume):.2f} µm³, {float(macro_pct):.2f}%")
    print(f"Mesopores (50–300 µm): {float(meso_volume):.2f} µm³, {float(meso_pct):.2f}%")
    print(f"Micropores (<50 µm): {float(micro_volume):.2f} µm³, {float(micro_pct):.2f}%")
    plt.bar(bin_centers, pore_size_freq, width=np.diff(bin_centers)[0], align='center')
    plt.title("Pore Size Distribution")
    plt.xlabel("Pore Diameter (µm)")
    plt.ylabel("Frequency")
    output_path = os.path.join(OUTPUT_DIR, "pore_size_distribution.png")
    plt.savefig(output_path)
    plt.close()
    if not os.path.exists(output_path):
        print(f"Error: Failed to save {output_path}")
    plt.hist(original_image[binary_image], bins=50, alpha=0.5, label='Pore Intensities', color='blue')
    plt.hist(original_image[~binary_image], bins=50, alpha=0.5, label='Solid Intensities', color='red')
    plt.axvline(FIXED_THRESHOLD, color='green', linestyle='--', label=f'Fixed Threshold = {FIXED_THRESHOLD}')
    plt.title("Thresholded Histogram")
    plt.xlabel("Intensity")
    plt.ylabel("Frequency")
    plt.legend()
    output_path = os.path.join(OUTPUT_DIR, "thresholded_histogram.png")
    plt.savefig(output_path)
    plt.close()
    if not os.path.exists(output_path):
        print(f"Error: Failed to save {output_path}")
    chi_density = compute_euler_number(binary_image)
    gamma = compute_gamma_indicator(binary_image)
    print(f"Euler Number Density (χ): {float(chi_density):.4f} mm⁻³")
    binary_image = morphology.remove_small_objects(binary_image, min_size=100)
    binary_image = morphology.remove_small_holes(binary_image, area_threshold=100)
    (num_pores, num_throats, avg_pore_diameter_net, avg_tortuosity, throat_size_freq, throat_bin_centers,
     avg_throat_diameter, throat_radii, tortuosity_values) = analyze_pore_network(binary_image)
    return (pore_volume, total_pore_area, max_diameter, min_diameter, avg_pore_sphericity,
            pore_surface_density, mean_curvature_density, percolates, critical_diameter,
            pore_size_freq, bin_centers, avg_diameter, pore_sizes, pore_diameters, pore_shape_factors,
            macro_volume, meso_volume, micro_volume, num_pores, num_throats, avg_pore_diameter_net,
            avg_tortuosity, throat_size_freq, throat_bin_centers, avg_throat_diameter, throat_radii, tortuosity_values,
            chi_density, gamma)

def compute_metrics(image, output_dir=None, voxel_size=None, fixed_threshold=None):
    """Compute pore analysis metrics, starting with fixed thresholding.

    Parameters are optional to preserve the original script behavior while avoiding hard-coded paths.
    """
    global OUTPUT_DIR, VOXEL_SIZE, VOXEL_VOLUME, FIXED_THRESHOLD
    if output_dir is not None:
        OUTPUT_DIR = str(output_dir)
        os.makedirs(OUTPUT_DIR, exist_ok=True)
    if voxel_size is not None:
        VOXEL_SIZE = np.array([voxel_size, voxel_size, voxel_size], dtype=np.float64) if np.isscalar(voxel_size) else np.array(voxel_size, dtype=np.float64)
        VOXEL_VOLUME = np.prod(VOXEL_SIZE)
    if fixed_threshold is not None:
        FIXED_THRESHOLD = fixed_threshold
    all_metrics = []
    print("\nProcessing Single Image")
    # Apply fixed threshold to segment pores (intensity < FIXED_THRESHOLD) from solids
    binary_image = image < FIXED_THRESHOLD
    if not np.any(binary_image):
        print(f"Warning: No pore voxels detected with fixed threshold {FIXED_THRESHOLD}. Check intensity range.")
    (pore_volume, total_pore_area, max_diameter, min_diameter, avg_pore_sphericity,
     pore_surface_density, mean_curvature_density, percolates, critical_diameter,
     pore_size_freq, bin_centers, avg_diameter, pore_sizes, pore_diameters, pore_shape_factors,
     macro_volume, meso_volume, micro_volume, num_pores, num_throats, avg_pore_diameter_net,
     avg_tortuosity, throat_size_freq, throat_bin_centers, avg_throat_diameter, throat_radii,
     tortuosity_values, chi_density, gamma) = analyze_pore_structure(binary_image, image)
    phi_vis = np.sum(binary_image, dtype=np.int64) / binary_image.size * 100
    metrics = {
        'Pore Volume (µm³)': float(pore_volume),
        'Total Pore Area (mm²)': float(total_pore_area),
        'Max Pore Diameter (µm)': float(max_diameter),
        'Min Pore Diameter (µm)': float(min_diameter),
        'Avg Pore Diameter (µm)': float(avg_diameter),
        'Avg Pore Sphericity': float(avg_pore_sphericity),
        'Pore Surface Density (mm⁻¹)': float(pore_surface_density),
        'Mean Curvature Density (mm⁻²)': float(mean_curvature_density),
        'Percolation': percolates,
        'Critical Pore Diameter (µm)': float(critical_diameter),
        'Macro Volume (µm³)': float(macro_volume),
        'Meso Volume (µm³)': float(meso_volume),
        'Micro Volume (µm³)': float(micro_volume),
        'Chi Density (mm⁻³)': float(chi_density),
        'Gamma': float(gamma),
        'Num Pores': num_pores,
        'Num Throats': num_throats,
        'Avg Pore Diameter (Network, µm)': float(avg_pore_diameter_net),
        'Avg Coordination Number': float(avg_tortuosity),  # Renamed from Avg Tortuosity
        'Avg Throat Diameter (µm)': float(avg_throat_diameter),
        'Phi_vis (%)': float(phi_vis),
        'Fixed Threshold': float(FIXED_THRESHOLD)
    }
    all_metrics.append(metrics)
    print(f"Debug: Saving pore_diameters.csv with {len(pore_sizes)} entries")
    output_path = os.path.join(OUTPUT_DIR, "pore_diameters.csv")
    try:
        # Save pore_diameters in chunks to handle large arrays
        chunk_size = 1000000
        with open(output_path, 'w') as f:
            f.write('Pore Diameter (µm)\n')
            for i in range(0, len(pore_sizes), chunk_size):
                chunk = pore_sizes[i:i + chunk_size]
                np.savetxt(f, chunk, fmt='%.6f')
        if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            print(f"Saved {output_path}")
        else:
            print(f"Error: Failed to save {output_path} or file is empty")
    except Exception as e:
        print(f"Error saving {output_path}: {e}")
    output_path = os.path.join(OUTPUT_DIR, "pore_size_distribution.csv")
    try:
        pd.DataFrame({
            'Pore Diameter (µm)': bin_centers,
            'Frequency': pore_size_freq
        }).to_csv(output_path, index=False)
        if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            print(f"Saved {output_path}")
        else:
            print(f"Error: Failed to save {output_path} or file is empty")
    except Exception as e:
        print(f"Error saving {output_path}: {e}")
    metrics_df = pd.DataFrame(all_metrics)
    print("\nMetrics for Single Image:")
    for key, value in metrics.items():
        if isinstance(value, float):
            print(f"{key}: {float(value):.6f}")
        else:
            print(f"{key}: {value}")
    output_path = os.path.join(OUTPUT_DIR, "analysis_results.xlsx")
    try:
        with pd.ExcelWriter(output_path) as writer:
            metrics_df.to_excel(writer, sheet_name="Metrics", index=False)
        if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            print(f"Saved {output_path}")
        else:
            print(f"Error: Failed to save {output_path} or file is empty")
    except Exception as e:
        print(f"Error saving {output_path}: {e}")
    return metrics_df


def run_pore_metrics(image_path, output_dir, voxel_size=15, fixed_threshold=85):
    """Run the complete pore-metrics workflow from an image path."""
    filtered_image = load_and_filter_image(image_path)
    return compute_metrics(filtered_image, output_dir=output_dir, voxel_size=voxel_size, fixed_threshold=fixed_threshold)
