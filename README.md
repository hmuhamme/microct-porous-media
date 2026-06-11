# microct-porous-media

[![PyPI version](https://img.shields.io/pypi/v/microct-porous-media.svg)](https://pypi.org/project/microct-porous-media/)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20649266.svg)](https://doi.org/10.5281/zenodo.20649266)

# microct-porous-media
A Python package for X-ray Microcomputed Tomography (µCT) image analysis of soils and horticultural substrates.

---

## Authors

**Hadi Hamaaziz Muhammed**
Faculty of Agricultural Sciences and Landscape Architecture
Osnabrück University of Applied Sciences, Germany
Email: [hadi.azizm@gmail.com](mailto:hadiazizm@gmail.com)

**Prof. Dr. Ruediger Anlauf**
Faculty of Agricultural Sciences and Landscape Architecture
Osnabrück University of Applied Sciences, Germany
Email: [r.anlauf@hs-osnabrueck.de](mailto:r.anlauf@hs-osnabrueck.de)

---

## Description

`microct-porous-media` is an open-source Python package for quantitative analysis of X-ray microcomputed tomography (µCT) images of soils and horticultural growing media.

The package integrates image preprocessing, representative elementary volume (REV) analysis, pore morphology characterization, pore network extraction, pore size distribution analysis, water retention curve estimation, permeability determination, and saturated hydraulic conductivity calculation into a unified workflow.

The package was developed from six original scientific µCT image-analysis scripts and reorganized into reusable, documented Python modules.

---

## Features

* Representative Elementary Volume (REV) analysis
* Pore morphology characterization
* Pore network extraction
* Critical pore diameter determination
* Pore size distribution analysis
* Water retention curve estimation
* Permeability calculation
* Saturated hydraulic conductivity estimation
* CSV and Excel export of results
* Publication-ready figures and plots

---

## Modules

### `rev_xy.py`

Two-dimensional Representative Elementary Volume (REV) analysis in the x-y plane, particularly suitable for horticultural substrates.

### `rev_xyz.py`

Three-dimensional REV analysis for mineral soils and full 3D characterization.

### `pore_metrics.py`

Computes:

* Total pore volume
* Pore size distribution
* Surface density
* Sphericity
* Euler number density (χ)
* Gamma indicator
* Percolation status
* Network metrics

Results can be exported to CSV and Excel formats.

### `critical_pore.py`

Determines critical pore and throat diameters using:

* PoreSpy network extraction
* Graph-based percolation-path analysis

### `pore_size_distribution.py`

Computes:

* Pore size distribution
* Macropore, mesopore, and micropore fractions
* Simplified water retention curve estimation using local thickness analysis

### `permeability.py`

Computes:

* Absolute permeability
* Saturated hydraulic conductivity

using:

* PoreSpy
* OpenPNM
* Stokes flow simulations

---

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

---

## Installation

### Install from PyPI

```bash
pip install microct-porous-media
```

### Install from Source

```bash
git clone https://github.com/hmuhamme/microct-porous-media.git

cd microct-porous-media

pip install -e .
```

### Conda Environment

```bash
conda env create -f environment.yml

conda activate microctpm

pip install -e .
```

---

## Documentation

A detailed user manual is included in the repository:

```text
docs/User_Manual_microct_porous_media_v1.0.pdf
```

and

```text
docs/User_Manual_microct_porous_media_v1.0.docx
```

The manual contains:

* Installation instructions
* Input requirements
* Parameter descriptions
* Example workflows
* Output interpretation
* Troubleshooting guidance

---

## Example

```python
from microctpm.permeability import calculate_permeability

result = calculate_permeability(
    image_path="sample.tif",
    output_dir="outputs/permeability",
    voxel_size=15e-6,
    threshold=85
)

print(result)
```

Example output:

```python
{
    'porosity': 0.618,
    'permeability_mD': 31981.79,
    'hydraulic_conductivity_cm_s': 0.0310
}
```

---

## Repository Structure

```text
microct-porous-media/

├── docs/
├── examples/
├── src/
│   └── microctpm/
├── tests/
├── README.md
├── LICENSE
├── INSTALL.md
├── requirements.txt
├── environment.yml
└── pyproject.toml
```

---

## Citation

If you use this package in research, please cite:

Muhammed, H.H., Anlauf, R., Daum, D., Schmitt Rahner, M., and Kuka, K. (2025).

**Pore Structure Analysis of Growing Media Using X-Ray Microcomputed Tomography.**

Applied Sciences, 15(22), 11886.

https://doi.org/10.3390/app152211886

---

## BibTeX

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

---

## Relationship to the Associated Publication

This package implements and extends the image-analysis workflow presented in the publication above, including:

* Representative Elementary Volume (REV) analysis
* Pore morphology characterization
* Pore network extraction
* Critical pore diameter estimation
* Water retention curve determination
* Permeability simulation
* Saturated hydraulic conductivity simulation using OpenPNM

Users are encouraged to cite the publication when using the methods implemented in this package.

---

## License

This project is distributed under the MIT License.

See:

```text
LICENSE
```

for details.

---
## Software DOI

This software archive is permanently available through Zenodo:

DOI: https://doi.org/10.5281/zenodo.20649266

Citation:

Muhammed, H. H., & Anlauf, R. (2026).
microct-porous-media (Version 1.0.2).
Zenodo.
https://doi.org/10.5281/zenodo.20649266

## Acknowledgements

The development of this package was supported by research activities at the Faculty of Agricultural Sciences and Landscape Architecture, Osnabrück University of Applied Sciences, Germany.

---

## Important Note

This package preserves the original scientific workflows as closely as possible while removing hard-coded paths and automatic execution.

Users are encouraged to validate results against experimental measurements and original research workflows before applying the package in scientific or engineering studies.

