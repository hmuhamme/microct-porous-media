# Migrated from original critical_pore_diameter.py
import numpy as np
import scipy.ndimage as ndi
import tifffile as tiff
import porespy as ps
import networkx as nx
import matplotlib.pyplot as plt
import scipy.sparse.csgraph as csgraph
import psutil

def trim_disconnected_pores(network):
    conns = network['throat.conns']
    pore_coords = network['pore.coords']
    num_pores = len(pore_coords)
    adj_matrix = np.zeros((num_pores, num_pores))
    adj_matrix[conns[:, 0], conns[:, 1]] = 1
    adj_matrix[conns[:, 1], conns[:, 0]] = 1
    n_components, labels = csgraph.connected_components(csgraph.csgraph_from_dense(adj_matrix, null_value=0), directed=False)
    if n_components > 1:
        print(f"Before trimming: {n_components} connected components")
        component_sizes = np.bincount(labels)
        largest_component = np.argmax(component_sizes)
        print(f"Largest component: {largest_component}, size: {component_sizes[largest_component]}")
        pores_to_keep = labels == largest_component
        throats_to_keep = np.all(pores_to_keep[conns], axis=1)
        network['pore.coords'] = pore_coords[pores_to_keep]
        network['throat.conns'] = network['throat.conns'][throats_to_keep]
        if 'throat.equivalent_diameter' in network:
            network['throat.equivalent_diameter'] = network['throat.equivalent_diameter'][throats_to_keep]
        if 'pore.equivalent_diameter' in network:
            network['pore.equivalent_diameter'] = network['pore.equivalent_diameter'][pores_to_keep]
        new_pore_indices = np.zeros(num_pores, dtype=int)
        new_pore_indices[pores_to_keep] = np.arange(np.sum(pores_to_keep))
        network['throat.conns'] = new_pore_indices[network['throat.conns']]
        print(f"After trimming: Number of pores = {len(network['pore.coords'])}, Number of throats = {len(network['throat.conns'])}")
    return network

def calculate_critical_pore_diameter(image_path, voxel_size=15e-6, threshold=110, r_max=10, sigma=1, output_dir=None, save_plot=True, show_plot=True):
    # Check available memory
    print(f"Available memory: {psutil.virtual_memory().available / 1e9:.2f} GB")

    # Load and preprocess image
    image = tiff.imread(image_path)
    print(f"Image shape: {image.shape}")

    fil_image = ndi.gaussian_filter(image, sigma=1)
    im = (fil_image < threshold).astype(np.uint8)
    im = ndi.binary_closing(im, structure=np.ones((3, 3, 3)))  # Stronger closing for connectivity

    # Pad image
    pad_z = (0, 1 if im.shape[0] % 2 != 0 else 0)
    pad_y = (0, 1 if im.shape[1] % 2 != 0 else 0)
    pad_x = (0, 1 if im.shape[2] % 2 != 0 else 0)
    im_padded = np.pad(im, pad_width=(pad_z, pad_y, pad_x), mode='constant', constant_values=1)
    print(f"Padded image shape: {im_padded.shape}")

    # Extract pore network
    try:
        print(f"Using r_max: {r_max} voxels ({r_max * voxel_size * 1e6:.2f} µm)")
        snow = ps.networks.snow2(im_padded, boundary_width=5, r_max=r_max, sigma=sigma)
        net = snow.network
        print(f"Shape after snow2: {snow.regions.shape}")
    except Exception as e:
        print(f"snow2 failed: {e}. Try reducing r_max or cropping the image.")
        return None

    # Trim disconnected pores
    net = trim_disconnected_pores(net)

    # Scale diameters
    if 'throat.equivalent_diameter' not in net:
        raise ValueError("No throat diameters available.")
    net['throat.equivalent_diameter'] *= voxel_size
    throat_diameters = net['throat.equivalent_diameter'] * 1e6  # Convert to µm
    throat_conns = net['throat.conns']
    pore_coords = net['pore.coords'] * voxel_size
    print(f"Throat diameters (µm): min={np.min(throat_diameters):.2f}, max={np.max(throat_diameters):.2f}")

    # Calculate critical pore diameter
    G = nx.Graph()
    max_diameter = np.max(throat_diameters) + 1
    for i, (p1, p2) in enumerate(throat_conns):
        G.add_edge(p1, p2, weight=max_diameter - throat_diameters[i])
    z_coords = pore_coords[:, 2]
    tol = voxel_size * 2
    inlet = np.where(z_coords < np.min(z_coords) + tol)[0]
    outlet = np.where(z_coords > np.max(z_coords) - tol)[0]
    print(f"Inlet pores: {len(inlet)}, Outlet pores: {len(outlet)}")

    min_diameters = []
    for start in inlet:
        for end in outlet:
            try:
                paths = nx.all_shortest_paths(G, start, end, weight='weight')
                for path in paths:
                    path_throats = []
                    for i in range(len(path) - 1):
                        p1, p2 = path[i], path[i + 1]
                        throat_idx = np.where((throat_conns[:, 0] == p1) & (throat_conns[:, 1] == p2) |
                                              (throat_conns[:, 0] == p2) & (throat_conns[:, 1] == p1))[0]
                        if len(throat_idx) > 0:
                            path_throats.append(throat_diameters[throat_idx[0]])
                    if path_throats:
                        min_diameters.append(np.min(path_throats))
            except nx.NetworkXNoPath:
                continue
            except nx.NodeNotFound:
                continue

    if not min_diameters:
        raise ValueError("No percolating path found from top to bottom of the sample.")

    critical_pore_diameter = max(min_diameters)
    print(f"Critical Pore Diameter: {critical_pore_diameter:.2f} µm")

    # Plot histogram
    plt.hist(throat_diameters, bins=50, edgecolor='k', alpha=0.7)
    plt.axvline(critical_pore_diameter, color='r', linestyle='--', label='Critical Diameter')
    plt.xlabel("Throat Diameter (µm)")
    plt.ylabel("Frequency")
    plt.legend()
    if save_plot:
        import os
        output_dir = os.getcwd() if output_dir is None else str(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        plt.savefig(os.path.join(output_dir, 'Throat_diameter_histogram.tif'), dpi=400)
    if show_plot:
        plt.show()
    else:
        plt.close()

    return critical_pore_diameter

# Run for the image
if __name__ == "__main__":
    try:
         calculate_critical_pore_diameter(
              "WF4S_1_cpd.tif",
            voxel_size=15e-6,
            threshold=110,
            r_max=10,
            sigma=1,
         )
    except MemoryError:
        print("MemoryError: Consider reducing r_max or cropping the image.")