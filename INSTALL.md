# Installation

```bash
cd microctpm_full_package
conda env create -f environment.yml
conda activate microctpm
pip install -e .
python -c "import microctpm; print(microctpm.__version__)"
```
