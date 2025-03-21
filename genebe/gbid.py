try:
    # try with a fast c-implementation ...
    import mmh3 as mmh3
except ImportError:
    # ... otherwise fallback to this code!
    import pymmh3 as mmh3


class InvalidPositionException(Exception):
    def __init__(self, message):
        super().__init__(message)


class PositionEncoder:
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
        return self.encode_1_based_without_throw(self, chromosome, position)

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
    BASES = ["A", "C", "G", "T"]

    positionEncoder = PositionEncoder()

    def hash_function(self, input_string: str, seed: int = 104729) -> int:
        # Compute the 128-bit hash from the input string using MurmurHash3 (hash128 returns a single 128-bit integer)
        hash_value = mmh3.hash128(input_string, seed)
        # Extract the high 64 bits and low 64 bits from the 128-bit hash
        low_part = hash_value & 0xFFFFFFFFFFFFFFFF  # Mask to get the lower 64 bits
        return low_part

    def set_value(self, input, mask, value):
        shift_count = (mask & -mask).bit_length() - 1
        new_value = (input & ~mask) | ((value << shift_count) & mask)
        return new_value

    def encode_bases(self, alt):
        bit_offset = 0
        encoded = 0
        for c in alt:
            value = 0
            if c == "A":
                value = 0
            elif c == "C":
                value = 1
            elif c == "G":
                value = 2
            elif c == "T":
                value = 3
            else:
                raise InvalidPositionException(
                    "Cannot encode alt as it contains strange letters: " + alt
                )

            encoded = self.set_value(encoded, 0b11 << bit_offset, value)

            bit_offset += 2

        return int(encoded)

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

    def encode1based(self, chromosome, position, ref, alt):
        # make 0 based
        position = position - 1

        altNorm = alt.upper()

        while len(ref) > 0 and len(altNorm) > 0 and ref[0] == altNorm[0]:
            # strip leading common character, may happen in deletions / insertions
            position += 1
            ref = ref[1:]
            altNorm = altNorm[1:]

        # if ref is string than replace with length
        if isinstance(ref, str):
            ref = len(ref)

        refLength = ref

        chromosomeNormalized = chromosome.lstrip("chr").upper()
        positionId = self.positionEncoder.encode_1_based_without_throw(
            chromosomeNormalized, position
        )

        if positionId == -1:
            raise InvalidPositionException(
                "Wrong position: {}:{}.".format(chromosome, position)
            )

        if not set(altNorm).issubset("ACGNT"):
            raise InvalidPositionException(
                "Alt should contain only ACGNT, position: {}, altNorm: {}.".format(
                    position, altNorm
                )
            )

        useChangeHash = False

        if len(altNorm) > self.MAX_INS_LENGTH or refLength > self.MAX_DEL_LENGTH:
            useChangeHash = True

        encodedAlt = 0

        if "N" in alt:
            useChangeHash = True
        else:
            encodedAlt = self.encode_bases(altNorm)

        if positionId == -1:
            return hash(chromosomeNormalized, position, refLength, altNorm)
        else:
            if useChangeHash:
                encoded = 0
                encoded = self.set_value(encoded, self.MASK_HASH, 0)
                encoded = self.set_value(encoded, self.MASK_C_HASH, 1)
                encoded = self.set_value(encoded, self.MASK_POS, positionId)

                changeHash = self.change_hash(refLength, altNorm)
                encoded = self.set_value(
                    encoded, self.MASK_CHANGE_HASH_VALUE, changeHash
                )

                return encoded

            else:
                encoded = 0
                # print('1'+bin(encoded))
                encoded = self.set_value(encoded, self.MASK_HASH, 0)
                # print('2'+bin(encoded))
                encoded = self.set_value(encoded, self.MASK_C_HASH, 0)
                # print('3'+bin(encoded))
                encoded = self.set_value(encoded, self.MASK_POS, positionId)
                # print('4'+bin(encoded))
                encoded = self.set_value(encoded, self.MASK_INS_LEN, len(altNorm))
                # print('5'+bin(encoded))
                encoded = self.set_value(encoded, self.MASK_DEL, refLength)
                # print('6'+bin(encoded))
                encoded = self.set_value(encoded, self.MASK_INS, encodedAlt)
                # print('7'+bin(encoded))
                return encoded


# Singleton instance of PositionEncoder
_variant_encoder_instance = VariantIdEncoder()


def encode_vcf_variant_gbid(chromosome, position, ref, alt):
    return _variant_encoder_instance.encode1based(chromosome, position, ref, alt)


# v = VariantIdEncoder()
# print(v.encode1based("1", 12, "A", "C"))
