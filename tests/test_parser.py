import pytest

parse_chunk_columns = pytest.importorskip("gtfread._parser").parse_chunk_columns


def test_known_and_dynamic_columns_are_parsed():
    lines = [
        (
            'gene_id "GENE1"; transcript_id "TX1"; gene_name "ABC1"; '
            'gene_type "protein_coding"; transcript_name "ABC1-201"; '
            'transcript_type "protein_coding"; exon_number "1"; exon_id "EXON1"; '
            'tag "basic"; havana_transcript "OTTT0001"; havana_gene "OTTG0001"; '
            'transcript_support_level "1"; hgnc_id "HGNC:5"; ccdsid "CCDS1"; '
            'artif_dupl "false"; level "2"; ont "PGO:0000001"; custom_key "custom";'
        ),
        'gene_id "GENE2"; gene_name "ABC1"; gene_type "protein_coding"; level "2";',
    ]

    columns = parse_chunk_columns(lines)

    assert columns["gene_id"] == ["GENE1", "GENE2"]
    assert columns["transcript_id"] == ["TX1", None]
    assert columns["havana_gene"] == ["OTTG0001", None]
    assert columns["hgnc_id"] == ["HGNC:5", None]
    assert columns["custom_key"] == ["custom", None]
    assert columns["gene_name"][0] is columns["gene_name"][1]
    assert columns["gene_type"][0] is columns["gene_type"][1]
    assert columns["level"][0] is columns["level"][1]
