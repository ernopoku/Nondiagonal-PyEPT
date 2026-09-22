# Installation


Use Python 3.10 or newer in a virtual environment:

```bash
git clone https://github.com/ernopoku/Nondiagonal-PyEPT.git
cd Nondiagonal-PyEPT
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
python examples/water.py
python -m nondiagonal_ept examples/hf.json -o hf_results.json
pytest -q
```

PySCF must have a wheel or a supported build on the chosen platform. Tested package versions are recorded in `docs/environment.txt`.

