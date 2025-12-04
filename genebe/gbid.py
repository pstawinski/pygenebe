class InvalidPositionException(Exception):
    def __init__(self, message):
        super().__init__(message)


class PositionEncoder:
    __slots__ = ("chromosomes", "chromosome_sizes", "sizes", "offsets")

    CHR_NOT_SUPPORTED = -1
    ERROR_WRONG_CHR_POSITION = -2

    def __init__(self):
        self.chromosomes = [
            "1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "9",
            "10",
            "11",
            "12",
            "13",
            "14",
            "15",
            "16",
            "17",
            "18",
            "19",
            "20",
            "21",
            "22",
            "X",
            "Y",
            "M",
        ]

        self.chromosome_sizes = {
            "1": 248956422,
            "2": 242193529,
            "3": 198295559,
            "4": 190214555,
            "5": 181538259,
            "6": 170805979,
            "7": 159345973,
            "8": 145138636,
            "9": 138394717,
            "10": 133797422,
            "11": 135086622,
            "12": 133275309,
            "13": 114364328,
            "14": 107043718,
            "15": 101991189,
            "16": 90338345,
            "17": 83257441,
            "18": 80373285,
            "19": 58617616,
            "20": 64444167,
            "21": 46709983,
            "22": 50818468,
            "X": 156040895,
            "Y": 57227415,
            "M": 16569,
        }
        self.sizes = [
            self.chromosome_sizes.get(chr_name, 0) for chr_name in self.chromosomes
        ]
        self.offsets = [sum(self.sizes[:i]) for i in range(len(self.chromosomes))]

    def encode(self, chromosome: str, position: int):
        return self.encode_1_based_without_throw(chromosome, position)

    def encode_1_based_without_throw(self, chromosome, position):
        chr_name = self._chr(chromosome)
        try:
            chr_index = self.chromosomes.index(chr_name)
        except ValueError:
            # chromosome not supported
            return self.CHR_NOT_SUPPORTED
        else:
            size = self.sizes[chr_index]
            if position > size or position < 0:
                print(f"Wrong position for {chromosome}: {position}")
                return self.ERROR_WRONG_CHR_POSITION
            offset = self.offsets[chr_index]
            return position - 1 + offset

    @staticmethod
    def _chr(chromosome):
        # This method corresponds to ChromosomeNameNormalizer.chr(chromosome) in the Java code
        return chromosome  # You may need to implement your normalization logic here


# Singleton instance of PositionEncoder
_position_encoder_instance = PositionEncoder()


def encode_vcf_position_gbid(chromosome, position):
    return _position_encoder_instance.encode_1_based_without_throw(chromosome, position)


class VariantIdEncoder:
    __slots__ = ("positionEncoder",)

    RADIX = 36
    MAX_DEL_LENGTH = (2**7) - 1
    MAX_INS_LENGTH = 9

    MASK_HASH = 0b1000000000000000000000000000000000000000000000000000000000000000
    MASK_NULL = 0b0100000000000000000000000000000000000000000000000000000000000000
    MASK_POS = 0b0011111111111111111111111111111111000000000000000000000000000000
    MASK_C_HASH = 0b0000000000000000000000000000000000100000000000000000000000000000
    MASK_INS_LEN = 0b0000000000000000000000000000000000011110000000000000000000000000
    MASK_DEL = 0b0000000000000000000000000000000000000001111111000000000000000000
    MASK_INS = 0b0000000000000000000000000000000000000000000000111111111111111111
    MASK_CHANGE_HASH_VALUE = (
        0b0000000000000000000000000000000000011111111111111111111111111111
    )

    # Pre-computed values for optimization
    BASE_ENCODING = {"A": 0, "C": 1, "G": 2, "T": 3}
    VALID_BASES = frozenset("ACGNT")

    # Pre-computed shift amounts
    SHIFT_POS = 30
    SHIFT_INS_LEN = 25
    SHIFT_DEL = 18
    SHIFT_C_HASH = 29

    def __init__(self):
        self.positionEncoder = PositionEncoder()

    def hash_function(self, input_string: str, seed: int = 104729) -> int:
        try:
            # try with a fast c-implementation ...
            import mmh3 as mmh3
        except ImportError:
            # ... otherwise fallback to this code!
            import pymmh3 as mmh3
        # Compute the 128-bit hash from the input string using MurmurHash3
        hash_value = mmh3.hash128(input_string, seed)
        # Extract the low 64 bits
        low_part = hash_value & 0xFFFFFFFFFFFFFFFF
        return low_part

    def set_value(self, input, mask, value):
        shift_count = (mask & -mask).bit_length() - 1
        new_value = (input & ~mask) | ((value << shift_count) & mask)
        return new_value

    def encode_bases_fast(self, alt):
        """Optimized base encoding using dict lookup and direct bit operations"""
        encoded = 0
        for i, c in enumerate(alt):
            value = self.BASE_ENCODING.get(c)
            if value is None:
                raise InvalidPositionException(
                    "Cannot encode alt as it contains strange letters: " + alt
                )
            encoded |= value << (i * 2)
        return encoded

    def hash(self, chromosome, position, del_length, alt):
        change = f"{chromosome}:{position}:{del_length}:{alt}"
        hash_value = self.hash_function(change.encode("utf-8"))
        value = self.set_value(hash_value, self.MASK_HASH, 1)
        return value

    def change_hash(self, del_length, alt):
        change = f"{del_length}:{alt}"
        hash_value = self.hash_function(change.encode("utf-8"))
        value = hash_value & self.MASK_CHANGE_HASH_VALUE
        return value

    def encodeSpdi(self, chromosome: str, position: int, del_length: int, ins: str):
        """
        Encode variant using SPDI format (0-based position, deletion length, insertion string)

        Args:
            chromosome: Chromosome name
            position: 0-based position
            del_length: Length of deletion (refLength)
            ins: Insertion string (normalized alt)
        """
        # Optimized chromosome normalization
        if chromosome.startswith(("chr", "CHR", "Chr")):
            chromosomeNormalized = chromosome[3:].upper()
        else:
            chromosomeNormalized = chromosome.upper()

        positionId = self.positionEncoder.encode_1_based_without_throw(
            chromosomeNormalized, position
        )

        if positionId == -1:
            raise InvalidPositionException(
                "Wrong position: {}:{}.".format(chromosome, position)
            )

        # Fast validation using frozenset
        if not self.VALID_BASES.issuperset(ins):
            raise InvalidPositionException(
                "Alt should contain only ACGNT, position: {}, altNorm: {}.".format(
                    position, ins
                )
            )

        # Check conditions once
        ins_len = len(ins)
        useChangeHash = (
            ins_len > self.MAX_INS_LENGTH or del_length > self.MAX_DEL_LENGTH
        )
        hasN = "N" in ins

        if useChangeHash or hasN:
            # Change hash path
            encoded = 0
            encoded = self.set_value(encoded, self.MASK_C_HASH, 1)
            encoded = self.set_value(encoded, self.MASK_POS, positionId)
            changeHash = self.change_hash(del_length, ins)
            encoded = self.set_value(encoded, self.MASK_CHANGE_HASH_VALUE, changeHash)
            return encoded
        else:
            # Direct encoding path - optimized with direct bit operations
            encodedAlt = self.encode_bases_fast(ins)

            # Combine all values in one operation
            encoded = (
                (positionId << self.SHIFT_POS)
                | (ins_len << self.SHIFT_INS_LEN)
                | (del_length << self.SHIFT_DEL)
                | encodedAlt
            )
            return encoded

    def encode1based(self, chromosome, position, ref, alt):
        # make 0 based
        position = position - 1

        altNorm = alt.upper()

        # Strip leading common characters
        while ref and altNorm and ref[0] == altNorm[0]:
            position += 1
            ref = ref[1:]
            altNorm = altNorm[1:]

        # if ref is string than replace with length
        refLength = len(ref) if isinstance(ref, str) else ref

        # Now call encodeSpdi with the SPDI format parameters
        return self.encodeSpdi(chromosome, position, refLength, altNorm)


# Singleton instance of VariantIdEncoder
_variant_encoder_instance = VariantIdEncoder()


def encode_vcf_variant_gbid(chromosome, position, ref, alt):
    """Optimized function to encode VCF variant"""
    return _variant_encoder_instance.encode1based(chromosome, position, ref, alt)


def encode_spdi_variant_gbid(chromosome, position0, del_length, ins):
    """Optimized function to encode SPDI variant"""
    return _variant_encoder_instance.encodeSpdi(chromosome, position0, del_length, ins)


# Example usage:
# v = VariantIdEncoder()
# print(v.encode1based("1", 12, "A", "C"))
