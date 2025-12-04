from enum import Enum
from functools import total_ordering


class InvalidGbIdentifierException(Exception):
    pass


class InvalidPositionException(Exception):
    pass


class GbIdentifierData:
    def __init__(self):
        self.type = None
        self.version = None
        self.transcript_id = None

    def set_type(self, type_val):
        self.type = type_val

    def set_version(self, version):
        self.version = version

    def set_transcript_id(self, transcript_id):
        self.transcript_id = transcript_id

    def __str__(self):
        return f"GbIdentifierData(type={self.type}, version={self.version}, transcript_id={self.transcript_id})"


class GbIdentifierType(Enum):
    ENST = ("ENST", 11, False, False)
    ENSP = ("ENSP", 11, False, False)
    ENSG = ("ENSG", 11, False, False)
    NM = ("NM_", 9, False, True)
    NR = ("NR_", 9, False, True)
    NP = ("NP_", 9, False, True)
    XM = ("XM_", 9, True, True)
    XP = ("XP_", 9, True, True)
    YP = ("YP_", 9, True, True)
    XR = ("XR_", 9, True, True)
    unassigned_transcript = ("unassigned_transcript_", 0, False, True)
    ENSR = ("ENSR", 11, False, False)

    def __init__(self, prefix, padding_size, computed, refseq):
        self.prefix = prefix
        self.padding_size = padding_size
        self.computed = computed
        self.refseq = refseq

    def get_prefix(self):
        return self.prefix

    def get_padding_size(self):
        return self.padding_size

    def is_computed(self):
        return self.computed

    def is_refseq(self):
        return self.refseq

    @classmethod
    def values(cls):
        return list(cls)


@total_ordering
class GbIdentifier:
    def __init__(self, data):
        if isinstance(data, str):
            self.data = GbIdentifierService.encode(data)
        else:
            self.data = data

    @staticmethod
    def parse(data):
        try:
            if isinstance(data, str):
                id_val = GbIdentifierService.encode(data)
                return GbIdentifier(id_val)
            else:
                return GbIdentifier(data)
        except Exception:
            return None

    def get_identifier(self):
        return GbIdentifierService.identifier(self.data)

    def get_version(self):
        return GbIdentifierService.version(self.data)

    def get_type(self):
        return GbIdentifierService.type(self.data)

    def get_data(self):
        return self.data

    def __str__(self):
        type_val = self.get_type()

        if type_val.is_refseq():
            # More complicated padding for RefSeq
            if self.get_type == GbIdentifierType.unassigned_transcript:
                padding_size = 0
            elif self.get_identifier() >= 1_000_000:
                padding_size = 9
            else:
                padding_size = 6
        else:
            padding_size = type_val.get_padding_size()

        result = type_val.get_prefix()
        result += str(self.get_identifier()).zfill(padding_size)

        version = self.get_version()
        if version > 0:
            result += f".{version}"

        return result

    def __hash__(self):
        return hash(self.data)

    def __eq__(self, other):
        if not isinstance(other, GbIdentifier):
            return False
        return self.data == other.data

    def __lt__(self, other):
        if not isinstance(other, GbIdentifier):
            return NotImplemented
        return self.data < other.data

    def equals_ignore_version(self, other):
        if not isinstance(other, GbIdentifier):
            return False
        return (
            self.get_identifier() == other.get_identifier()
            and self.get_type() == other.get_type()
        )


class GbIdentifierService:
    # Bit masks for different parts of the identifier
    MASK_OMIT = 0b1000000000000000000000000000000000000000000000000000000000000000
    MASK_TYPE = 0b0111111110000000000000000000000000000000000000000000000000000000
    MASK_ID = 0b0000000001111111111111111111111111111111111111000000000000000000
    MASK_VERSION = 0b0000000000000000000000000000000000000000000000111111111111111111

    @staticmethod
    def decode(id_val):
        pid = GbIdentifierData()
        type_id = GbIdentifierService._read_value(id_val, GbIdentifierService.MASK_TYPE)
        pid.set_type(GbIdentifierService._type_id_to_string(type_id))
        pid.set_version(
            GbIdentifierService._read_value(id_val, GbIdentifierService.MASK_VERSION)
        )
        pid.set_transcript_id(
            GbIdentifierService._read_value_long(id_val, GbIdentifierService.MASK_ID)
        )
        return pid

    @staticmethod
    def without_version(id_val):
        return id_val & ~GbIdentifierService.MASK_VERSION

    @staticmethod
    def version(id_val):
        return GbIdentifierService._read_value(id_val, GbIdentifierService.MASK_VERSION)

    @staticmethod
    def type(id_val):
        type_id = GbIdentifierService._read_value(id_val, GbIdentifierService.MASK_TYPE)
        return GbIdentifierService._type_id_to_string(type_id)

    @staticmethod
    def to_string(id_val):
        if id_val is None:
            return None
        else:
            return str(GbIdentifier.parse(id_val))

    @staticmethod
    def identifier(id_val):
        return GbIdentifierService._read_value_long(id_val, GbIdentifierService.MASK_ID)

    @staticmethod
    def equal(left, right):
        if isinstance(left, GbIdentifier):
            return left.get_data() == right
        return left == right

    @staticmethod
    def equal_without_version(left, right):
        if isinstance(left, GbIdentifier):
            left = left.get_data()
        return GbIdentifierService.without_version(
            left
        ) == GbIdentifierService.without_version(right)

    @staticmethod
    def encode(transcript):
        transcript_id = 0
        version = 0
        real_type = None

        # Find matching type by prefix
        for gb_type in GbIdentifierType.values():
            if transcript.startswith(gb_type.get_prefix()):
                real_type = gb_type
                break

        if real_type is None:
            raise InvalidGbIdentifierException(f"Wrong type for {transcript}")

        if "." in transcript:
            # Extract transcript ID between prefix and "."
            start_idx = len(real_type.get_prefix())
            dot_idx = transcript.find(".")
            transcript_id = int(transcript[start_idx:dot_idx])
            version = int(transcript[dot_idx + 1 :])
        else:
            # Extract transcript ID after prefix
            transcript_id = int(transcript[len(real_type.get_prefix()) :])

        type_id = list(GbIdentifierType).index(real_type)

        return GbIdentifierService._encode(transcript_id, version, type_id)

    @staticmethod
    def _encode(transcript_id, version, type_id):
        encoded = 0
        encoded = GbIdentifierService._set_value(
            encoded, GbIdentifierService.MASK_TYPE, type_id
        )
        encoded = GbIdentifierService._set_value_long(
            encoded, GbIdentifierService.MASK_ID, transcript_id
        )
        encoded = GbIdentifierService._set_value(
            encoded, GbIdentifierService.MASK_VERSION, version
        )
        return encoded

    @staticmethod
    def _type_id_to_string(id_val):
        return list(GbIdentifierType)[id_val]

    @staticmethod
    def _set_value(input_val, mask, value):
        shift_count = (mask & -mask).bit_length() - 1  # Count trailing zeros
        new_value = (input_val & ~mask) | (((value << shift_count) & mask))
        return new_value

    @staticmethod
    def _set_value_long(input_val, mask, value):
        shift_count = (mask & -mask).bit_length() - 1  # Count trailing zeros
        new_value = (input_val & ~mask) | (((value << shift_count) & mask))
        return new_value

    @staticmethod
    def _read_value(input_val, mask):
        shift_count = (mask & -mask).bit_length() - 1  # Count trailing zeros
        return int((input_val & mask) >> shift_count)

    @staticmethod
    def _read_value_long(input_val, mask):
        shift_count = (mask & -mask).bit_length() - 1  # Count trailing zeros
        return (input_val & mask) >> shift_count


def encode_transcript_id(transcript):
    """
    Easy-to-use wrapper for encoding transcript ID string to long.

    Args:
        transcript (str): Transcript ID string (e.g., "ENST00000404276.6", "NM_001234567.3")

    Returns:
        int: Encoded long value

    Raises:
        InvalidGbIdentifierException: If transcript format is invalid
    """
    return GbIdentifierService.encode(transcript)


def decode_transcript_id(encoded_id):
    """
    Easy-to-use wrapper for decoding long value back to transcript ID string.

    Args:
        encoded_id (int): Encoded long value

    Returns:
        str: Transcript ID string (e.g., "ENST00000404276.6", "NM_001234567.3")
    """
    return GbIdentifierService.to_string(encoded_id)


# Example usage
if __name__ == "__main__":
    try:
        # Test the new easy-to-use wrapper functions
        print("=== Testing wrapper functions ===")

        # Test ENST transcript
        transcript1 = "ENST00000404276.6"
        encoded1 = encode_transcript_id(transcript1)
        decoded1 = decode_transcript_id(encoded1)
        print(f"Original: {transcript1}")
        print(f"Encoded:  {encoded1}")
        print(f"Decoded:  {decoded1}")
        print(f"Round trip successful: {transcript1 == decoded1}")

        # Test RefSeq transcript
        transcript2 = "NM_001234567.3"
        encoded2 = encode_transcript_id(transcript2)
        decoded2 = decode_transcript_id(encoded2)
        print(f"\nOriginal: {transcript2}")
        print(f"Encoded:  {encoded2}")
        print(f"Decoded:  {decoded2}")
        print(f"Round trip successful: {transcript2 == decoded2}")

        print("\n=== Testing original functionality ===")

        # Test encoding and decoding
        id_val = GbIdentifierService.encode("ENST00000404276.6")
        print(f"Encoded long value: {id_val}")
        decoded = GbIdentifierService.decode(id_val)
        print(f"Decoded: {decoded}")

        # Test with GbIdentifier from string
        gb_id = GbIdentifier.parse("ENST00000404276.6")
        print(f"\nGbIdentifier from string: {gb_id}")
        print(f"Internal long: {gb_id.get_data()}")
        print(f"Type: {gb_id.get_type()}")
        print(f"Identifier: {gb_id.get_identifier()}")
        print(f"Version: {gb_id.get_version()}")

        # Test creating GbIdentifier from the long value (round trip)
        gb_id_from_long = GbIdentifier.parse(gb_id.get_data())
        print(f"\nGbIdentifier from long {gb_id.get_data()}: {gb_id_from_long}")
        print(f"Round trip successful: {gb_id == gb_id_from_long}")

        # Test other identifier types
        refseq_id = GbIdentifier.parse("NM_001234567.3")
        print(f"\nRefSeq ID: {refseq_id}")
        print(f"RefSeq internal long: {refseq_id.get_data()}")
        print(f"RefSeq Type: {refseq_id.get_type()}")

        # Test RefSeq with large identifier (should use 9-digit padding)
        large_refseq = GbIdentifier.parse("NM_1234567890.1")
        print(f"\nLarge RefSeq: {large_refseq}")
        print(f"Large RefSeq internal long: {large_refseq.get_data()}")

        # Test round trip with large RefSeq
        large_refseq_from_long = GbIdentifier.parse(large_refseq.get_data())
        print(f"Large RefSeq from long: {large_refseq_from_long}")
        print(f"Large RefSeq round trip: {large_refseq == large_refseq_from_long}")

    except Exception as e:
        print(f"Error: {e}")
