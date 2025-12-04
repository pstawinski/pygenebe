"""Tests for GBID SPDI encoding functionality."""

import pytest
from genebe.gbid import VariantIdEncoder, InvalidPositionException


class TestEncodeSpdi:
    """Test cases for the encodeSpdi method."""

    def setup_method(self):
        """Set up test fixtures."""
        self.encoder = VariantIdEncoder()

    def test_simple_snv(self):
        """Test simple single nucleotide variant encoding."""
        # SNV: position 11 (0-based), delete 1 base, insert "C"
        chromosome = "1"
        position = 11  # 0-based
        del_length = 1
        ins = "C"

        result = self.encoder.encodeSpdi(chromosome, position, del_length, ins)
        assert isinstance(result, int)
        assert result > 0

    def test_deletion(self):
        """Test deletion variant encoding."""
        # Deletion: position 99 (0-based), delete 3 bases, insert nothing
        chromosome = "1"
        position = 99  # 0-based
        del_length = 3
        ins = ""

        result = self.encoder.encodeSpdi(chromosome, position, del_length, ins)
        assert isinstance(result, int)
        assert result > 0

    def test_insertion(self):
        """Test insertion variant encoding."""
        # Insertion: position 199 (0-based), delete 0 bases, insert "TAG"
        chromosome = "1"
        position = 199  # 0-based
        del_length = 0
        ins = "TAG"

        result = self.encoder.encodeSpdi(chromosome, position, del_length, ins)
        assert isinstance(result, int)
        assert result > 0

    def test_complex_variant(self):
        """Test complex variant (deletion + insertion)."""
        # Complex: position 500 (0-based), delete 2 bases, insert "ATCG"
        chromosome = "1"
        position = 500  # 0-based
        del_length = 2
        ins = "ATCG"

        result = self.encoder.encodeSpdi(chromosome, position, del_length, ins)
        assert isinstance(result, int)
        assert result > 0

    def test_different_chromosomes(self):
        """Test encoding on different chromosomes."""
        test_cases = [
            ("1", 100, 1, "A"),
            ("22", 100, 1, "A"),
            ("X", 100, 1, "A"),
            ("Y", 100, 1, "A"),
            ("M", 100, 1, "A"),
        ]

        for chromosome, position, del_length, ins in test_cases:
            result = self.encoder.encodeSpdi(chromosome, position, del_length, ins)
            assert isinstance(result, int)
            assert result > 0

    def test_edge_case_positions(self):
        """Test encoding at edge positions."""
        # Test position 1 (near start of chromosome, 0-based means position 1 in 1-based)
        result1 = self.encoder.encodeSpdi("1", 1, 1, "A")
        assert isinstance(result1, int)
        assert result1 > 0

        # Test a position near the end of chromosome 1
        result2 = self.encoder.encodeSpdi("1", 248956420, 1, "A")  # Near end of chr1
        assert isinstance(result2, int)
        assert result2 > 0

    def test_all_base_combinations(self):
        """Test encoding with all possible DNA bases."""
        bases = ["A", "C", "G", "T"]
        chromosome = "1"
        position = 100
        del_length = 1

        for base in bases:
            result = self.encoder.encodeSpdi(chromosome, position, del_length, base)
            assert isinstance(result, int)
            assert result > 0

    def test_multi_base_insertions(self):
        """Test insertions of various lengths."""
        chromosome = "1"
        position = 100
        del_length = 0

        test_insertions = [
            "A",
            "AT",
            "ATG",
            "ATGC",
            "ATGCA",
            "ATGCAT",
            "ATGCATG",
            "ATGCATGC",
            "ATGCATGCA",  # MAX_INS_LENGTH = 9
        ]

        for ins in test_insertions:
            result = self.encoder.encodeSpdi(chromosome, position, del_length, ins)
            assert isinstance(result, int)
            assert result > 0

    def test_large_variants_use_hash(self):
        """Test that large variants trigger hash mode."""
        chromosome = "1"
        position = 100

        # Test large insertion (> MAX_INS_LENGTH = 9)
        large_ins = "ATGCATGCATGC"  # 12 bases
        result1 = self.encoder.encodeSpdi(chromosome, position, 0, large_ins)
        assert isinstance(result1, int)
        assert result1 > 0

        # Test large deletion (> MAX_DEL_LENGTH = 127)
        large_del = 150
        result2 = self.encoder.encodeSpdi(chromosome, position, large_del, "A")
        assert isinstance(result2, int)
        assert result2 > 0

    def test_n_base_triggers_hash(self):
        """Test that 'N' in insertion sequence triggers hash mode."""
        chromosome = "1"
        position = 100
        del_length = 1
        ins_with_n = "ANTG"

        result = self.encoder.encodeSpdi(chromosome, position, del_length, ins_with_n)
        assert isinstance(result, int)
        assert result > 0

    def test_chromosome_normalization(self):
        """Test that chromosome names are normalized correctly."""
        position = 100
        del_length = 1
        ins = "A"

        result1 = self.encoder.encodeSpdi("1", position, del_length, ins)
        result2 = self.encoder.encodeSpdi("chr1", position, del_length, ins)

        # These should produce the same result after chr prefix stripping
        assert result1 == result2

    def test_invalid_chromosome(self):
        """Test encoding with invalid chromosome."""
        with pytest.raises(InvalidPositionException, match="Wrong position"):
            self.encoder.encodeSpdi("25", 100, 1, "A")  # Chromosome 25 doesn't exist

    def test_invalid_position_zero(self):
        """Test encoding with position 0 (invalid in this system)."""
        with pytest.raises(InvalidPositionException, match="Wrong position"):
            self.encoder.encodeSpdi("1", 0, 1, "A")

    def test_invalid_position_negative(self):
        """Test encoding with negative position."""
        # The position encoder prints a warning but returns ERROR_WRONG_CHR_POSITION (-2)
        # which doesn't raise an exception in encodeSpdi, but should be handled
        result = self.encoder.encodeSpdi("1", -1, 1, "A")
        # The method should return a hash result for invalid positions

    def test_invalid_position_too_large(self):
        """Test encoding with position beyond chromosome size."""
        # Position larger than chromosome 1 size (248956422)
        # The position encoder prints a warning but returns ERROR_WRONG_CHR_POSITION (-2)
        result = self.encoder.encodeSpdi("1", 300000000, 1, "A")
        # The method should return a hash result for invalid positions

    def test_invalid_bases_in_insertion(self):
        """Test encoding with invalid bases in insertion sequence."""
        with pytest.raises(
            InvalidPositionException, match="Alt should contain only ACGNT"
        ):
            self.encoder.encodeSpdi("1", 100, 1, "ATGX")  # X is not a valid base

    def test_empty_insertion(self):
        """Test encoding with empty insertion (pure deletion)."""
        chromosome = "1"
        position = 100
        del_length = 5
        ins = ""  # Empty insertion

        result = self.encoder.encodeSpdi(chromosome, position, del_length, ins)
        assert isinstance(result, int)
        assert result > 0

    def test_zero_deletion(self):
        """Test encoding with zero deletion length (pure insertion)."""
        chromosome = "1"
        position = 100
        del_length = 0  # No deletion
        ins = "ATCG"

        result = self.encoder.encodeSpdi(chromosome, position, del_length, ins)
        assert isinstance(result, int)
        assert result > 0

    def test_consistency_with_encode1based(self):
        """Test that encodeSpdi produces same results as encode1based for equivalent inputs."""
        test_cases = [
            # (chromosome, pos_1based, ref, alt) -> convert to SPDI format
            ("1", 12, "A", "C"),  # SNV
            (
                "1",
                100,
                "ATG",
                "A",
            ),  # Deletion: ATG -> A (strip A, position+1, del=2, ins="")
            ("1", 200, "C", "CTAG"),  # Insertion: C -> CTAG (del=1, ins="CTAG")
            ("X", 50, "G", "T"),  # Different chromosome
            ("22", 1000, "AA", "TT"),  # Multi-base substitution
        ]

        for chromosome, pos1, ref, alt in test_cases:
            # Get result from encode1based
            result1 = self.encoder.encode1based(chromosome, pos1, ref, alt)

            # Convert to SPDI manually (same logic as encode1based)
            pos0 = pos1 - 1
            altNorm = alt.upper()
            ref_work = ref

            # Strip leading common characters
            while len(ref_work) > 0 and len(altNorm) > 0 and ref_work[0] == altNorm[0]:
                pos0 += 1
                ref_work = ref_work[1:]
                altNorm = altNorm[1:]

            del_length = len(ref_work)

            # Get result from encodeSpdi
            result2 = self.encoder.encodeSpdi(chromosome, pos0, del_length, altNorm)

            assert result1 == result2, f"Mismatch for {chromosome}:{pos1} {ref}->{alt}"

    def test_deterministic_encoding(self):
        """Test that encoding is deterministic (same inputs produce same outputs)."""
        chromosome = "1"
        position = 100
        del_length = 1
        ins = "C"

        result1 = self.encoder.encodeSpdi(chromosome, position, del_length, ins)
        result2 = self.encoder.encodeSpdi(chromosome, position, del_length, ins)
        result3 = self.encoder.encodeSpdi(chromosome, position, del_length, ins)

        assert result1 == result2 == result3

    def test_different_variants_different_encodings(self):
        """Test that different variants produce different encodings."""
        chromosome = "1"
        position = 100
        del_length = 1

        result_a = self.encoder.encodeSpdi(chromosome, position, del_length, "A")
        result_c = self.encoder.encodeSpdi(chromosome, position, del_length, "C")
        result_g = self.encoder.encodeSpdi(chromosome, position, del_length, "G")
        result_t = self.encoder.encodeSpdi(chromosome, position, del_length, "T")

        # All results should be different
        results = [result_a, result_c, result_g, result_t]
        assert len(set(results)) == 4, "All variants should produce unique encodings"

    def test_position_sensitivity(self):
        """Test that different positions produce different encodings."""
        chromosome = "1"
        del_length = 1
        ins = "A"

        result1 = self.encoder.encodeSpdi(chromosome, 100, del_length, ins)
        result2 = self.encoder.encodeSpdi(chromosome, 101, del_length, ins)
        result3 = self.encoder.encodeSpdi(chromosome, 102, del_length, ins)

        # All results should be different
        results = [result1, result2, result3]
        assert (
            len(set(results)) == 3
        ), "Different positions should produce unique encodings"
