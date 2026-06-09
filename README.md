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

A Python package for Microcomputed Tomography X-Ray (µCT) image analysis of soils and horticultural substrates, including representative elementary volume (REV) analysis, pore morphology characterization, pore network extraction, pore size distribution, water retention curve estimation, permeability, and saturated hydraulic conductivity determination.

A Python package migrated from six original µCT image-analysis research scripts for soils and horticultural substrates.

## Modules

- `rev_xy.py`: two-dimensional REV analysis in the x-y plane, recommended for horticultural substrates.
- `rev_xyz.py`: three-dimensional REV analysis, recommended for mineral soils or 3D REV assessment.
- `pore_metrics.py`: pore volume, pore size, surface density, sphericity, Euler number density (χ), gamma indicator, percolation, network metrics, and CSV/Excel outputs.
- `critical_pore.py`: critical pore/throat diameter using PoreSpy network extraction and graph-based percolating-path analysis.
- `pore_size_distribution.py`: pore size distribution and simplified water-retention curve from local thickness.
- `permeability.py`: permeability and saturated hydraulic conductivity using PoreSpy + OpenPNM + Stokes flow.

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

### Option 1: Python virtual environment

```bash
python -m venv venv
venv\Scripts\activate
pip install -U pip
pip install -e .
```

### Option 2: Conda environment

```bash
conda env create -f environment.yml
conda activate microctpm
pip install -e .
```



## Example

```python
from microctpm.permeability import calculate_permeability

result = calculate_permeability(
    image_path="PeatS_2.tif",
    output_dir="outputs/permeability",
    voxel_size=15e-6,
    threshold=85,
)
print(result)
```
## Citation

If you use this package in research, please cite:

Muhammed, H.H., Anlauf, R., Daum, D., Schmitt Rahner, M., and Kuka, K. (2025).

**Pore Structure Analysis of Growing Media Using X-Ray Microcomputed Tomography.**

Applied Sciences, 15(22), 11886.

https://doi.org/10.3390/app152211886

### BibTeX

```bibtex
@article{Muhammed2025MicroCT,
  author = {Muhammed, Hadi Hamaaziz and Anlauf, Ruediger and Daum, Diemo and Schmitt Rahner, Mayka and Kuka, Katrin},
  title = {Pore Structure Analysis of Growing Media Using X-Ray Microcomputed Tomography},
  journal = {Applied Sciences},
  year = {2025},
  volume = {15},
  number = {22},
  pages = {11886},
  doi = {10.3390/app152211886}
}
```

### Relationship to this package

This package implements and extends the image-analysis workflow presented in the publication above, including:

- Representative Elementary Volume (REV) analysis
- Pore morphology characterization
- Pore network extraction
- Critical pore diameter estimation
- Water retention curve (WRC) determination
- Saturated hydraulic conductivity simulation using OpenPNM

Users are encouraged to cite the publication when using the methods implemented in this package.
## Important note

This version preserves the original scientific workflows as closely as possible while removing hard-coded paths and automatic execution. Validate all results against the original scripts before public release.
