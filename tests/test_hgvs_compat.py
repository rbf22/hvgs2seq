"""Tests for the hgvs_compat module."""
import pytest
from unittest.mock import MagicMock, patch
from hgvs2seq.hgvs_compat import (
    SequenceVariant, SimplePosition, Interval, PosEdit,
    HGVSParseError, Parser, VariantMapper, exceptions
)

# Test data
SIMPLE_HGVS_G = "NC_000001.11:g.12345678G>A"
SIMPLE_HGVS_C = "NM_001.2:c.123A>G"
SIMPLE_HGVS_N = "NR_123.4:n.56+7A>G"
SIMPLE_HGVS_P = "NP_001.2:p.Val600Glu"

# Fixtures
@pytest.fixture
def simple_position():
    return SimplePosition(base=123, uncertain=False)

@pytest.fixture
def interval(simple_position):
    return Interval(start=simple_position, end=simple_position, uncertain=False)

@pytest.fixture
def posedit(interval):
    return PosEdit(pos=interval, edit=None)

@pytest.fixture
def genomic_variant(posedit):
    return SequenceVariant(ac="NC_000001.11", type="g", posedit=posedit)

@pytest.fixture
def coding_variant(posedit):
    return SequenceVariant(ac="NM_001.2", type="c", posedit=posedit)

# Tests
class TestSimplePosition:
    def test_initialization(self):
        pos = SimplePosition(base=123, uncertain=True)
        assert pos.base == 123
        assert pos.uncertain is True

class TestInterval:
    def test_initialization(self, simple_position):
        interval = Interval(start=simple_position, end=simple_position, uncertain=True)
        assert interval.start == simple_position
        assert interval.end == simple_position
        assert interval.uncertain is True

class TestPosEdit:
    def test_initialization(self, interval):
        posedit = PosEdit(pos=interval, edit="A>G")
        assert posedit.pos == interval
        assert posedit.edit == "A>G"

class TestSequenceVariant:
    def test_initialization(self, posedit):
        variant = SequenceVariant(ac="NM_001.2", type="c", posedit=posedit)
        assert variant.ac == "NM_001.2"
        assert variant.type == "c"
        assert variant.posedit == posedit

    def test_invalid_variant_type(self, posedit):
        with pytest.raises(ValueError, match="Invalid variant type: x"):
            SequenceVariant(ac="NM_001.2", type="x", posedit=posedit)

    def test_str_representation(self, posedit):
        variant = SequenceVariant(ac="NM_001.2", type="c", posedit=posedit)
        assert str(variant) == "NM_001.2:c.123_123"

class TestParser:
    @pytest.fixture
    def parser(self):
        return Parser()

    def test_parse_genomic_variant(self, parser):
        variant = parser.parse_hgvs_variant(SIMPLE_HGVS_G)
        assert variant.ac == "NC_000001.11"
        assert variant.type == "g"
        assert variant.posedit.pos.start.base == 12345678

    def test_parse_coding_variant(self, parser):
        variant = parser.parse_hgvs_variant(SIMPLE_HGVS_C)
        assert variant.ac.startswith("NM_")
        assert variant.type == "c"

    def test_parse_noncoding_variant(self, parser):
        variant = parser.parse_hgvs_variant(SIMPLE_HGVS_N)
        assert variant.ac.startswith("NR_")
        assert variant.type == "n"

    def test_parse_protein_variant(self, parser):
        variant = parser.parse_hgvs_variant(SIMPLE_HGVS_P)
        assert variant.ac.startswith("NP_")
        assert variant.type == "p"

    def test_parse_invalid_variant(self, parser):
        with pytest.raises(HGVSParseError):
            parser.parse_hgvs_variant("invalid_hgvs_string")

class TestVariantMapper:
    @pytest.fixture
    def mapper(self):
        return VariantMapper(hdp=None)

    def test_g_to_c_conversion(self, mapper, genomic_variant):
        coding_variant = mapper.g_to_c(genomic_variant, "NM_001.2")
        assert coding_variant.type == "c"
        assert coding_variant.ac == "NM_001.2"
        assert coding_variant.posedit == genomic_variant.posedit

    def test_c_to_g_conversion(self, mapper, coding_variant):
        genomic_variant = mapper.c_to_g(coding_variant, "NC_000001.11")
        assert genomic_variant.type == "g"
        assert genomic_variant.ac == "NC_000001.11"
        assert genomic_variant.posedit == coding_variant.posedit

    def test_n_to_c_conversion(self, mapper):
        noncoding_variant = SequenceVariant(
            ac="NR_123.4",
            type="n",
            posedit=PosEdit(
                pos=Interval(
                    start=SimplePosition(base=56, uncertain=False),
                    end=SimplePosition(base=56, uncertain=False),
                    uncertain=False
                ),
                edit=None
            )
        )
        coding_variant = mapper.n_to_c(noncoding_variant, "NM_001.2")
        assert coding_variant.type == "c"
        assert coding_variant.ac == "NM_001.2"

    def test_c_to_c_conversion(self, mapper, coding_variant):
        new_coding_variant = mapper.c_to_c(coding_variant, "NM_002.3")
        assert new_coding_variant.type == "c"
        assert new_coding_variant.ac == "NM_002.3"
        assert new_coding_variant.posedit == coding_variant.posedit

    def test_invalid_conversion(self, mapper, genomic_variant):
        # Test invalid conversion (g to n)
        with pytest.raises(ValueError, match="Expected n. variant, got g. variant"):
            mapper.n_to_c(genomic_variant, "NM_001.2")

class TestExceptions:
    def test_hgvs_parse_error(self):
        with pytest.raises(HGVSParseError):
            raise HGVSParseError("Test error")

    def test_hgvs_data_not_available_error(self):
        with pytest.raises(exceptions.HGVSDataNotAvailableError):
            raise exceptions.HGVSDataNotAvailableError("Test error")
