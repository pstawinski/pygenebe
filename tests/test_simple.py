import pytest
from genebe.transcriptid import GbIdentifier, GbIdentifierService, GbIdentifierType


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
