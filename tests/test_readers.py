from pathlib import Path

import pandas as pd
from pandas.testing import assert_frame_equal

from gtfread import (
    find_first_data_line_index,
    read_gtf,
    read_gtf_full,
    read_gtf_full_python,
    read_gtf_python,
    read_gtf_restricted,
    read_gtf_restricted_python,
)


def _normalize_frame_for_compare(df):
    df = df.copy()
    for column in ["Chromosome", "Source", "Feature", "Strand", "Frame"]:
        if column in df.columns:
            df[column] = df[column].astype(str)
    return df.astype(object).where(pd.notna(df), None)


def _write_temp_gtf(tmp_path: Path, contents: str) -> Path:
    path = tmp_path / "test.gtf"
    path.write_text(contents)
    return path


def test_find_first_data_line_index_skips_comments(tmp_path: Path):
    path = _write_temp_gtf(tmp_path, "# comment\n## comment\nchr1\tsrc\tgene\t1\t2\t.\t+\t.\tgene_id \"G1\";\n")
    assert find_first_data_line_index(path) == 2


def test_read_gtf_parses_semicolons_inside_quoted_attributes(tmp_path: Path):
    path = _write_temp_gtf(
        tmp_path,
        "chr1\tRefSeq\tCDS\t20486313\t20486315\t.\t+\t0\t"
        'gene_id "SELENOO"; transcript_id "NM_001115017.5"; db_xref "GeneID:417745"; '
        'gbkey "CDS"; gene "SELENOO"; note "UGA stop codon recoded as selenocysteine; '
        'The RefSeq protein has 1 substitution compared to this genomic sequence"; '
        'product "protein adenylyltransferase SelO, mitochondrial"; '
        'protein_id "NP_001108489.5"; transl_except "(pos:20486313..20486315,aa:Sec)"; '
        'exon_number "1";\n',
    )

    result = read_gtf(path)
    row = result.iloc[0]

    assert row["note"] == (
        "UGA stop codon recoded as selenocysteine; "
        "The RefSeq protein has 1 substitution compared to this genomic sequence"
    )
    assert row["product"] == "protein adenylyltransferase SelO, mitochondrial"
    assert row["exon_number"] == "1"
    assert row["Start"] == 20486312


def test_read_gtf_duplicate_attr_keeps_all_values(tmp_path: Path):
    path = _write_temp_gtf(
        tmp_path,
        "chr1\thavana\texon\t11869\t12227\t.\t+\t.\t"
        'gene_id "ENSG1"; transcript_id "ENST1"; exon_id "ENSE1"; tag "CCDS"; tag "basic";\n',
    )

    default = read_gtf(path)
    duplicate = read_gtf(path, duplicate_attr=True)

    assert default.iloc[0]["tag"] == "basic"
    assert duplicate.iloc[0]["tag"] == "CCDS,basic"


def test_read_gtf_full_multiple_chunks_match(tmp_path: Path):
    path = _write_temp_gtf(
        tmp_path,
        "# header\n"
        'chr1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG1"; gene_name "DDX11L1"; level "2";\n'
        'chr1\thavana\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; transcript_name "DDX11L1-201"; tag "basic";\n'
        'chr1\thavana\texon\t12010\t12057\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; exon_number "1"; exon_id "ENSE1";\n',
    )

    skiprows = find_first_data_line_index(path)
    chunked = read_gtf_full(path, skiprows=skiprows, chunksize=2)
    unchunked = read_gtf_full(path, skiprows=skiprows, chunksize=1000)

    chunked = _normalize_frame_for_compare(chunked)
    unchunked = _normalize_frame_for_compare(unchunked)

    assert_frame_equal(chunked, unchunked, check_dtype=False)


def test_read_gtf_restricted_returns_core_columns(tmp_path: Path):
    path = _write_temp_gtf(
        tmp_path,
        "# header\n"
        'chr1\thavana\texon\t12010\t12057\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; exon_number "1"; exon_id "ENSE1";\n',
    )

    skiprows = find_first_data_line_index(path)
    result = read_gtf_restricted(path, skiprows=skiprows)

    assert list(result.columns) == [
        "Chromosome",
        "Source",
        "Feature",
        "Start",
        "End",
        "Score",
        "Strand",
        "Frame",
        "gene_id",
        "transcript_id",
        "exon_number",
        "exon_id",
    ]
    assert result.iloc[0]["Start"] == 12009
    assert result.iloc[0]["transcript_id"] == "ENST1"


def test_read_gtf_python_matches_compiled_reader(tmp_path: Path):
    path = _write_temp_gtf(
        tmp_path,
        "# header\n"
        'chr1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id "ENSG1"; gene_name "DDX11L1"; level "2";\n'
        'chr1\thavana\ttranscript\t11869\t14409\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; transcript_name "DDX11L1-201"; tag "basic";\n',
    )

    skiprows = find_first_data_line_index(path)
    compiled = read_gtf_full(path, skiprows=skiprows, chunk_size=2)
    python = read_gtf_full_python(path, skiprows=skiprows, chunk_size=2)

    compiled = _normalize_frame_for_compare(compiled)
    python = _normalize_frame_for_compare(python)

    assert list(python["gene_id"]) == list(compiled["gene_id"])
    assert list(python["gene_name"]) == ["DDX11L1", None]
    assert list(python["transcript_id"]) == [None, "ENST1"]
    assert list(python["transcript_name"]) == [None, "DDX11L1-201"]
    assert list(python["tag"]) == [None, "basic"]
    assert list(python["Start"]) == list(compiled["Start"])


def test_read_gtf_python_duplicate_attr_keeps_all_values(tmp_path: Path):
    path = _write_temp_gtf(
        tmp_path,
        "chr1\thavana\texon\t11869\t12227\t.\t+\t.\t"
        'gene_id "ENSG1"; transcript_id "ENST1"; exon_id "ENSE1"; tag "CCDS"; tag "basic";\n',
    )

    compiled = read_gtf(path, duplicate_attr=True)
    python = read_gtf_python(path, duplicate_attr=True)

    assert compiled.iloc[0]["tag"] == "CCDS,basic"
    assert python.iloc[0]["tag"] == "CCDS,basic"


def test_read_gtf_restricted_python_matches_compiled(tmp_path: Path):
    path = _write_temp_gtf(
        tmp_path,
        "# header\n"
        'chr1\thavana\texon\t12010\t12057\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1"; exon_number "1"; exon_id "ENSE1";\n',
    )

    skiprows = find_first_data_line_index(path)
    compiled = read_gtf_restricted(path, skiprows=skiprows)
    python = read_gtf_restricted_python(path, skiprows=skiprows)

    compiled = _normalize_frame_for_compare(compiled)
    python = _normalize_frame_for_compare(python)

    assert_frame_equal(compiled, python, check_dtype=False)
