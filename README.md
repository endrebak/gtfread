# gtfread

`gtfread` is a small Cython-backed package for parsing and reading GTF files into pandas dataframes.

## Layout

- `src/gtfread/_parser.pyx`: the Cython attribute parser
- `src/gtfread/_fallback.py`: pure Python fallback for source-tree imports before the extension is built
- `src/gtfread/readers.py`: the dataframe reader API
- `tests/test_parser.py`: parser smoke test
- `tests/test_readers.py`: GTF reader tests

## Install

```bash
python -m pip install -e .
```

That uses `pyproject.toml` to install the build requirements and compile the Cython extension.

## Usage

```python
from gtfread import read_gtf

df = read_gtf("annotation.gtf")
```

## Build a wheel

```bash
python -m build
```

## In-place extension build

```bash
python setup.py build_ext --inplace
```

## Do you need `setup.py`?

No. A modern project can build a Cython extension from `pyproject.toml` alone if you declare the extension there and use a build backend that supports it. `setup.py` is still useful when you want `cythonize(...)` directly, custom build logic, or more control over compiler options.
