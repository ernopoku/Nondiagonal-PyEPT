# Installation

Use Python 3.10 or newer and a supported PySCF platform.

```bash
git clone https://github.com/ernopoku/Nondiagonal-PyEPT.git
cd Nondiagonal-PyEPT
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
pytest -q
python examples/water.py
python -m nondiagonal_ept examples/hf.json -o hf_results.json
```

Dependencies are NumPy, SciPy, and PySCF. The test extra installs pytest. On Windows, use the activation command appropriate to your environment and a supported PySCF installation.

[Recorded validation environment](https://github.com/ernopoku/Nondiagonal-PyEPT/blob/main/docs/environment.txt)

[[Home]] | [[Examples]]
