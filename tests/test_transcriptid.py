"""Tests for transcript ID encoding and decoding functionality using pytest."""

import pytest
from genebe.transcriptid import (
    GbIdentifier,
    GbIdentifierService,
    GbIdentifierType,
    InvalidGbIdentifierException,
)


def test_basic_encoding():
    """Test basic encoding functionality."""
    transcript_id = "ENST00000404276.6"
    expected_long = 105978527750

    actual_long = GbIdentifierService.encode(transcript_id)
    assert actual_long == expected_long


def test_gb_identifier_creation():
    """Test creating GbIdentifier objects."""
    transcript_id = "ENST00000404276.6"
    expected_long = 105978527750

    gb_id = GbIdentifier(transcript_id)
    assert gb_id.get_data() == expected_long


def test_round_trip():
    """Test that encoding and decoding works correctly."""
    transcript_id = "ENST00000404276.6"

    encoded_long = GbIdentifierService.encode(transcript_id)
    gb_id = GbIdentifier(encoded_long)
    decoded_string = str(gb_id)

    assert decoded_string == transcript_id


def test_equality():
    """Test equality comparison."""
    gb_id1 = GbIdentifier("ENST00000404276.6")
    gb_id2 = GbIdentifier("ENST00000404276.6")
    gb_id3 = GbIdentifier("ENST00000123456.1")

    assert gb_id1 == gb_id2
    assert gb_id1 != gb_id3


def test_service_methods():
    """Test static methods of GbIdentifierService."""
    test_long = 105978527750  # ENST00000404276.6

    version = GbIdentifierService.version(test_long)
    assert version == 6

    gb_type = GbIdentifierService.type(test_long)
    assert gb_type == GbIdentifierType.ENST

    identifier = GbIdentifierService.identifier(test_long)
    assert identifier == 404276

    string_repr = GbIdentifierService.to_string(test_long)
    assert string_repr == "ENST00000404276.6"


@pytest.fixture
def test_cases():
    """Test cases with expected transcript ID to long value mappings."""
    return [
        ("ENST00000404276.6", 105978527750),
        ("ENST00000123456.1", 32363249665),
        ("NM_001234567.3", 108086714691223555),
        ("NM_123456.2", 108086423420141570),
        ("ENST00000404276", 105978527744),
        ("ENST00000123456", 32363249664),
        ("NM_001234567", 108086714691223552),
        ("NM_123456", 108086423420141568),
    ]


def test_encode_transcript_ids(test_cases):
    """Test encoding of transcript identifiers to long values."""
    for transcript_id, expected_long in test_cases:
        actual_long = GbIdentifierService.encode(transcript_id)
        assert actual_long == expected_long


def test_round_trip_encoding_decoding(test_cases):
    """Test that encoding and then decoding returns the original string."""
    for transcript_id, _ in test_cases:
        encoded_long = GbIdentifierService.encode(transcript_id)
        gb_id = GbIdentifier(encoded_long)
        decoded_string = str(gb_id)
        assert decoded_string == transcript_id


def test_invalid_identifiers():
    """Test handling of invalid identifier formats."""
    with pytest.raises(InvalidGbIdentifierException):
        GbIdentifierService.encode("INVALID00000123456.1")

    with pytest.raises((InvalidGbIdentifierException, ValueError)):
        GbIdentifierService.encode("")


@pytest.mark.parametrize(
    "transcript_id,expected_type",
    [
        ("ENST00000123456.1", GbIdentifierType.ENST),
        ("ENSP00000123456.1", GbIdentifierType.ENSP),
        ("ENSG00000123456.1", GbIdentifierType.ENSG),
        ("NM_123456.1", GbIdentifierType.NM),
    ],
)
def test_gb_identifier_types(transcript_id, expected_type):
    """Test all supported GbIdentifierType enums."""
    gb_id = GbIdentifier(transcript_id)
    assert gb_id.get_type() == expected_type
