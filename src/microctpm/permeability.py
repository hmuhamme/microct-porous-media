
"""Permeability and saturated hydraulic conductivity from pore networks.

Migrated from the original script Permeability_Conductivity.py. The hard-coded image
path and immediate execution have been replaced by the callable function
``calculate_permeability``. The core workflow is preserved: image filtering,
thresholding, PoreSpy snow2 extraction, OpenPNM conversion, connected-component
trimming, geometry scaling, Valvatne-Blunt hydraulic conductance, Stokes flow,
absolute permeability, and saturated hydraulic conductivity.
"""
from __future__ import annotations

from pathlib import Path
import os
import numpy as np
import scipy.ndimage as ndi
import scipy.sparse.csgraph as csgraph
import tifffile as tiff
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    import porespy as ps
except Exception:  # pragma: no cover
    ps = None
try:
    import openpnm as op
except Exception:  # pragma: no cover
    op = None


def trim_disconnected_pores(network):
    """Trim an OpenPNM network to its largest connected pore component."""
    if op is None:
        raise ImportError("openpnm is required for trim_disconnected_pores")
    conns = network['throat.conns']
    adj_matrix = np.zeros((network.Np, network.Np))
    adj_matrix[conns[:, 0], conns[:, 1]] = 1
    adj_matrix[conns[:, 1], conns[:, 0]] = 1
    n_components, labels = csgraph.connected_components(
        csgraph.csgraph_from_dense(adj_matrix, null_value=0), directed=False
    )
    if n_components > 1:
        print(f"Before trimming: {n_components} connected components")
        component_sizes = np.bincount(labels)
        largest_component = np.argmax(component_sizes)
        print(f"Largest component: {largest_component}, size: {component_sizes[largest_component]}")
        pores_to_keep = labels == largest_component
        throats_to_keep = np.all(pores_to_keep[conns], axis=1)
        op.topotools.trim(network=network, pores=np.where(~pores_to_keep)[0], throats=np.where(~throats_to_keep)[0])
        print(f"After trimming: Number of pores = {network.Np}, Number of throats = {network.Nt}")
    return network


def calculate_permeability(
    image_path,
    output_dir=None,
    voxel_size=15e-6,
    threshold=85,
    sigma=1,
    boundary_width=5,
    inlet_pressure=2000.0,
    outlet_pressure=1990.0,
    viscosity=1.0e-3,
    density=1000.0,
    gravity=9.81,
    flow_axis=0,
    tol=None,
    save_outputs=True,
    make_plot=True,
    show_plot=True,
    plot_filename='3d_pore_network.tif',
    save_plots=None,
    show_plots=None,
    pressure_drop=None
):
    """Calculate absolute permeability and saturated hydraulic conductivity.

    Parameters
    ----------
    image_path : str or Path
        3D tif/tiff µCT image.
    output_dir : str or Path, optional
        Folder for outputs. If None, current working directory is used.
    voxel_size : float
        Voxel spacing in metres, e.g. 15e-6 or 46e-6.
    threshold : float
        Intensity threshold. Voxels below threshold are pores.
    sigma : float
        Gaussian filter sigma.
    boundary_width : int
        PoreSpy snow2 boundary width.
    inlet_pressure, outlet_pressure : float
        Pressure boundary conditions in Pa.
    flow_axis : int
        Axis along which flow is simulated. 0=z/depth in array indexing, 1=y, 2=x.

    Returns
    -------
    dict
        Permeability and hydraulic conductivity plus diagnostic values.
    """
    if ps is None:
        raise ImportError("porespy is required. Install with conda or pip before running permeability.")
    if op is None:
        raise ImportError("openpnm is required. Install with conda or pip before running permeability.")

    # Backward-compatible aliases used by some test scripts
    if save_plots is not None:
        make_plot = bool(save_plots)
        save_outputs = bool(save_plots)
    if show_plots is not None:
        show_plot = bool(show_plots)

    output_dir = Path.cwd() if output_dir is None else Path(output_dir)
    if save_outputs:
        output_dir.mkdir(parents=True, exist_ok=True)

    image = tiff.imread(str(image_path))
    fil_image = ndi.gaussian_filter(image, sigma=sigma)
    print(f"Image intensity range: min={np.min(fil_image):.2f}, max={np.max(fil_image):.2f}")
    print(f"Image shape: {image.shape}")
    im = (fil_image < threshold).astype(np.uint8)
    im = ndi.binary_closing(im, structure=np.ones((1, 1, 1)))
    porosity = np.mean(im)
    print(f"Porosity: {porosity:.3f}")

    if make_plot and save_outputs:
        plt.figure()
        plt.hist(fil_image.ravel(), bins=256, alpha=0.5, label='Image')
        plt.axvline(x=threshold, color='g', linestyle='--', label=f'Threshold={threshold:.2f}')
        plt.title('Intensity Histogram')
        plt.legend()
        plt.savefig(output_dir / 'intensity_histogram.tif', dpi=400)
        if show_plot:
            plt.show()
        else:
            plt.close()
        plt.figure()
        plt.imshow(im[im.shape[0]//2, :, :], cmap='gray')
        plt.title(f'Binary Image Slice (Depth={im.shape[0]//2}, thresh={threshold:.2f})')
        plt.savefig(output_dir / 'binary_image_slice.tif', dpi=400)
        if show_plot:
            plt.show()
        else:
            plt.close()

    pad_z = (0, 1 if im.shape[0] % 2 != 0 else 0)
    pad_y = (0, 1 if im.shape[1] % 2 != 0 else 0)
    pad_x = (0, 1 if im.shape[2] % 2 != 0 else 0)
    im_padded = np.pad(im, pad_width=(pad_z, pad_y, pad_x), mode='constant', constant_values=1)
    print(f"Padded image shape for snow2: {im_padded.shape}")

    snow = ps.networks.snow2(im_padded, boundary_width=boundary_width)
    net = snow.network
    print(f"Shape after snow2 processing: {snow.regions.shape}")

    pn = op.io.network_from_porespy(net)
    pn = trim_disconnected_pores(pn)

    spacing = float(voxel_size)
    dims = snow.regions.shape
    pn['pore.coords'] *= spacing
    if 'pore.equivalent_diameter' in pn:
        pn['pore.equivalent_diameter'] *= spacing
    if 'throat.equivalent_diameter' in pn:
        pn['throat.equivalent_diameter'] *= spacing

    pn.add_model_collection(op.models.collections.geometry.spheres_and_cylinders)
    pn.regenerate_models()

    pn['pore.shape_factor'] = np.ones(pn.Np) >= 0.048
    pn['throat.shape_factor'] = np.ones(pn.Nt) >= 0.048
    pore_radius = pn['pore.equivalent_diameter'] / 2
    throat_radius = pn['throat.equivalent_diameter'] / 2
    pn['pore.area'] = np.pi * pore_radius**2
    pn['throat.cross_sectional_area'] = np.pi * throat_radius**2
    throat_length = pn['throat.length']
    L1 = throat_length * 0.5
    L2 = throat_length * 0.5
    Lt = throat_length
    pn['throat.conduit_lengths'] = np.vstack((L1, Lt, L2)).T

    print("Pore diameter stats (µm): min =", np.min(pn['pore.equivalent_diameter'])*1e6,
          ", max =", np.max(pn['pore.equivalent_diameter'])*1e6,
          ", mean =", np.mean(pn['pore.equivalent_diameter'])*1e6)
    print("Throat diameter stats (µm): min =", np.min(pn['throat.equivalent_diameter'])*1e6,
          ", max =", np.max(pn['throat.equivalent_diameter'])*1e6,
          ", mean =", np.mean(pn['throat.equivalent_diameter'])*1e6)

    coords = pn['pore.coords']
    axis_coords = coords[:, flow_axis]
    if tol is None:
        tol = 2 * spacing
    inlet = pn.pores()[axis_coords < (axis_coords.min() + tol)]
    outlet = pn.pores()[axis_coords > (axis_coords.max() - tol)]
    print(f"Number of inlet pores: {len(inlet)}")
    print(f"Number of outlet pores: {len(outlet)}")
    if len(inlet) == 0 or len(outlet) == 0:
        raise ValueError("Inlet or outlet pores are empty. Check network extraction or adjust tolerance.")

    phase = op.phase.Phase(network=pn, name='water')
    phase['pore.viscosity'] = viscosity
    phase['throat.viscosity'] = viscosity
    phase.add_model_collection(op.models.collections.physics.basic)
    phase.add_model(propname='throat.hydraulic_conductance',
                    model=op.models.physics.hydraulic_conductance.valvatne_blunt,
                    pore_viscosity='pore.viscosity',
                    throat_viscosity='throat.viscosity',
                    pore_shape_factor='pore.shape_factor',
                    throat_shape_factor='throat.shape_factor',
                    pore_area='pore.area',
                    throat_area='throat.cross_sectional_area',
                    conduit_lengths='throat.conduit_lengths')
    phase.regenerate_models()

    flow = op.algorithms.StokesFlow(network=pn, phase=phase)
    flow.set_value_BC(pores=inlet, values=inlet_pressure)
    flow.set_value_BC(pores=outlet, values=outlet_pressure)
    flow.run()
    phase.update(flow.soln)

    Q = flow.rate(pores=inlet, mode='group')[0]
    L = dims[flow_axis] * spacing
    other_axes = [i for i in range(3) if i != flow_axis]
    A = (dims[other_axes[0]] * spacing) * (dims[other_axes[1]] * spacing)
    if pressure_drop is None:
        delta_p = abs(inlet_pressure - outlet_pressure)
    else:
        delta_p = float(pressure_drop)
    K = Q * L * viscosity / (A * delta_p)

    
    K_mD = K / 0.987e-12 * 1000
    k_cm_s = K * (density * gravity / viscosity) * 100
    print(f'The permeability is: {K_mD:.2f} mD')
    print(f'The hydraulic conductivity is: {k_cm_s:.4f} cm/s')

    if make_plot:
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection="3d")

        op.visualization.plot_connections(pn, ax=ax)

        relative_pressure = phase["pore.pressure"] - np.min(phase["pore.pressure"])

        coords_plot = pn["pore.coords"]
        pressure_scatter = ax.scatter(
          coords_plot[:, 0],
          coords_plot[:, 1],
          coords_plot[:, 2],
          c=relative_pressure,
          s=8,
          cmap="viridis",
    )

    ax.set_xlabel("X (m)", fontsize=12, labelpad=15)
    ax.set_ylabel("Y (m)", fontsize= 12, labelpad=15)
    ax.set_zlabel("Z (m)", fontsize= 12, labelpad=20)
    ax.tick_params(axis='x', pad=5, labelsize=10)
    ax.tick_params(axis='y', pad=5, labelsize=10)
    ax.tick_params(axis='z', pad=8, labelsize=10)

    fig.colorbar(
        pressure_scatter,
        ax=ax,
        label="Relative pressure (Pa)",
        pad=0.10,
        shrink=0.6,
    )

    if save_outputs:
        plt.tight_layout()
        plt.savefig(output_dir / plot_filename, dpi=400, bbox_inches="tight")

    if show_plot:
        plt.show()
    else:
        plt.close(fig)

    return {
        'porosity': float(porosity),
        'flow_rate_m3_s': float(Q),
        'permeability_m2': float(K),
        'permeability_mD': float(K_mD),
        'hydraulic_conductivity_cm_s': float(k_cm_s),
        'num_pores': int(pn.Np),
        'num_throats': int(pn.Nt),
        'num_inlet_pores': int(len(inlet)),
        'num_outlet_pores': int(len(outlet)),
        'pressure_drop_Pa': float(delta_p)
    }


def calculate_hydraulic_conductivity(*args, **kwargs):
    """Alias for calculate_permeability."""
    return calculate_permeability(*args, **kwargs)
