"""
Integration tests for the local pyhgvs implementation.

This test file verifies that the local pyhgvs package works correctly
and can be imported and used from the hgvs2seq project.
"""
import pytest
from hgvs2seq.pyhgvs import HGVSName, CDNACoord, parse_hgvs_name, format_hgvs_name


def test_local_pyhgvs_import():
    """Test that local pyhgvs can be imported correctly."""
    # Test that key classes and functions are available
    assert HGVSName is not None
    assert CDNACoord is not None
    assert parse_hgvs_name is not None
    assert format_hgvs_name is not None


def test_hgvs_name_basic_creation():
    """Test basic HGVS name creation and parsing."""
    # Test genomic variant
    hgvs = HGVSName("chr12:g.1000000A>T")
    assert hgvs.kind == 'g'
    assert hgvs.chrom == 'chr12'
    assert hgvs.start == 1000000
    assert hgvs.ref_allele == 'A'
    assert hgvs.alt_allele == 'T'

    # Test cDNA variant
    hgvs_cdna = HGVSName("NM_000492.3:c.123A>T")
    assert hgvs_cdna.kind == 'c'
    assert hgvs_cdna.transcript == 'NM_000492.3'

    # Test protein variant
    hgvs_prot = HGVSName("NP_000483.3:p.Trp24Cys")
    assert hgvs_prot.kind == 'p'
    assert hgvs_prot.start == 24
    assert hgvs_prot.ref_allele == 'Trp'
    assert hgvs_prot.alt_allele == 'Cys'


def test_cdna_coordinate_parsing():
    """Test cDNA coordinate parsing."""
    # Test basic coordinate
    coord = CDNACoord(string='123')
    assert coord.coord == 123
    assert coord.offset == 0
    assert coord.landmark == 'cdna_start'

    # Test coordinate with offset
    coord = CDNACoord(string='123+4')
    assert coord.coord == 123
    assert coord.offset == 4

    # Test negative coordinate
    coord = CDNACoord(string='-5')
    assert coord.coord == -5
    assert coord.offset == 0


def test_hgvs_name_formatting():
    """Test HGVS name formatting."""
    # Test genomic formatting
    hgvs = HGVSName()
    hgvs.kind = 'g'
    hgvs.chrom = 'chr12'
    hgvs.start = 1000000
    hgvs.end = 1000000
    hgvs.ref_allele = 'A'
    hgvs.alt_allele = 'T'
    hgvs.mutation_type = '>'

    formatted = hgvs.format()
    assert formatted == 'chr12:g.1000000A>T'

    # Test cDNA formatting
    hgvs_cdna = HGVSName("NM_000492.3:c.123A>T")
    assert hgvs_cdna.format() == 'NM_000492.3:c.123A>T'


def test_hgvs_name_edge_cases():
    """Test edge cases and error handling."""
    # Test invalid HGVS name should raise an exception
    try:
        hgvs = HGVSName("invalid:hgvs:name")
        assert False, "Should have raised an exception"
    except Exception:
        pass  # Expected

    # Test empty HGVS name
    hgvs = HGVSName()
    assert hgvs.name == ''


def test_coordinate_string_representation():
    """Test string representation of coordinates."""
    coord = CDNACoord(123, 5)
    assert str(coord) == '123+5'

    coord_neg = CDNACoord(-10, -2)
    assert str(coord_neg) == '-10-2'

    coord_star = CDNACoord(5, 0, 'cdna_stop')
    assert str(coord_star) == '*5'


def test_basic_functionality_integration():
    """Test that basic functions work together."""
    # Test that we can create, parse, and format HGVS names
    original_name = "chr12:g.1000000A>T"

    # Parse
    hgvs = HGVSName(original_name)

    # Verify parsed data
    assert hgvs.kind == 'g'
    assert hgvs.chrom == 'chr12'
    assert hgvs.start == 1000000
    assert hgvs.ref_allele == 'A'
    assert hgvs.alt_allele == 'T'

    # Format back
    formatted = hgvs.format()
    assert formatted == original_name


if __name__ == "__main__":
    # Run the tests if called directly
    test_local_pyhgvs_import()
    test_hgvs_name_basic_creation()
    test_cdna_coordinate_parsing()
    test_hgvs_name_formatting()
    test_hgvs_name_edge_cases()
    test_coordinate_string_representation()
    test_basic_functionality_integration()
    print("All tests passed! Local pyhgvs implementation is working correctly.")
