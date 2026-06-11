# microct-porous-media

## Authors

**Hadi Hamaaziz Muhammed**  
Faculty of Agricultural Sciences and Landscape Architecture  
Osnabrück University of Applied Sciences  
Email: hadiazizm@gmail.com

**Prof. Dr. Ruediger Anlauf**  
Faculty of Agricultural Sciences and Landscape Architecture  
Osnabrück University of Applied Sciences  
Email: r.anlauf@hs-osnabrueck.de

## Description

`microct-porous-media` is a Python package for X-ray microcomputed tomography (µCT) image analysis of soils and horticultural substrates.

The package supports:

- Representative elementary volume (REV) analysis
- Pore morphology characterization
- Pore network extraction
- Pore size distribution analysis
- Water retention curve estimation
- Permeability calculation
- Saturated hydraulic conductivity determination

The package was migrated from six original µCT image-analysis research scripts and reorganized into reusable Python modules.

## Modules

- `rev_xy.py`: two-dimensional REV analysis in the x-y plane, recommended for horticultural substrates.
- `rev_xyz.py`: three-dimensional REV analysis, recommended for mineral soils or full 3D REV assessment.
- `pore_metrics.py`: pore volume, pore size, surface density, sphericity, Euler number density (χ), gamma indicator, percolation, network metrics, and CSV/Excel outputs.
- `critical_pore.py`: critical pore/throat diameter using PoreSpy network extraction and graph-based percolating-path analysis.
- `pore_size_distribution.py`: pore size distribution and simplified water-retention curve from local thickness.
- `permeability.py`: permeability and saturated hydraulic conductivity using PoreSpy, OpenPNM, and Stokes flow.

## Scientific Workflow

The package provides a complete workflow for µCT image analysis:

1. Image preprocessing and filtering
2. REV analysis (2D and 3D)
3. Pore morphology characterization
4. Critical pore diameter determination
5. Pore size distribution analysis
6. Water retention curve estimation
7. Pore network extraction
8. Permeability determination
9. Saturated hydraulic conductivity calculation

## Installation

### Install from PyPI

```bash
pip install microct-porous-media
