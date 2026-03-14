# gtfread

`gtfread` is a small package for parsing and reading GTF files into pandas dataframes.

## Layout

- `gtfread/_parser.pyx`: the optional compiled attribute parser
- `gtfread/readers.py`: the serial dataframe reader API
- `tests/test_parser.py`: parser smoke test
- `tests/test_readers.py`: GTF reader tests

## Install

```bash
python -m pip install -e .
```

That uses `pyproject.toml` to install the build requirements and compile the Cython extension.

## Usage

```python
from gtfread import read_gtf, read_gtf_python

df = read_gtf("annotation.gtf")
df_python = read_gtf_python("annotation.gtf")
```

`read_gtf(...)` uses the compiled parser path when available. `read_gtf_python(...)` uses the high-level pure Python parser path used by the current `pyrunges` reader style.

If you want to use the compiled low-level parser directly, pass it raw attribute strings from column 9 of the GTF before they have been expanded:

```python
import pandas as pd

from gtfread import find_first_data_line_index, parse_chunk_columns

skiprows = find_first_data_line_index("annotation.gtf")
attribute_lines = pd.read_csv(
    "annotation.gtf",
    sep="\t",
    header=None,
    usecols=[8],
    names=["Attribute"],
    skiprows=skiprows,
)["Attribute"].tolist()

compiled_columns = parse_chunk_columns(attribute_lines)
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

No. A modern project can build a Cython extension from `pyproject.toml` alone if you declare the extension there and use a build backend that supports it. `setup.py` is still useful when you want `cythonize(...)` directly, custom build logic, or more control over compiler options. If the compiled extension is unavailable, the pure Python reader still works.
