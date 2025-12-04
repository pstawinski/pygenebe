#!/usr/bin/env python3
"""
Demonstration of the new encodeSpdi method for GBID encoding.

This shows the difference between encode1based (VCF format) and encodeSpdi (SPDI format),
and demonstrates various use cases for the new SPDI-based encoding method.
"""

from genebe.gbid import VariantIdEncoder


def main():
    print("=== GBID encodeSpdi Demonstration ===\n")

    encoder = VariantIdEncoder()

    print("1. Basic comparison: encode1based vs encodeSpdi")
    print("-" * 50)

    # Example 1: Simple SNV
    chr1, pos1, ref, alt = "1", 12, "A", "C"
    result1 = encoder.encode1based(chr1, pos1, ref, alt)

    # Equivalent SPDI format: 0-based position, deletion length, insertion
    pos0 = pos1 - 1  # Convert to 0-based
    del_length = len(ref)
    ins = alt
    result2 = encoder.encodeSpdi(chr1, pos0, del_length, ins)

    print(f"VCF format:  {chr1}:{pos1} {ref}->{alt}")
    print(f"encode1based result: {result1}")
    print(f"SPDI format: {chr1}:{pos0} del={del_length} ins='{ins}'")
    print(f"encodeSpdi result:   {result2}")
    print(f"Results match: {result1 == result2}")
    print()

    print("2. Complex variant examples")
    print("-" * 30)

    examples = [
        # (description, chr, pos_1based, ref, alt)
        ("Simple SNV", "1", 100, "G", "T"),
        ("Deletion", "2", 200, "ATG", "A"),
        ("Insertion", "3", 300, "C", "CTAG"),
        ("Complex variant", "X", 400, "GCA", "TTAG"),
        ("Large insertion", "1", 500, "A", "ATCGATCGATCG"),  # Triggers hash mode
    ]

    for desc, chrom, pos1, ref, alt in examples:
        print(f"{desc}: {chrom}:{pos1} {ref}->{alt}")

        # encode1based
        result1 = encoder.encode1based(chrom, pos1, ref, alt)

        # Convert to SPDI format manually (same logic as in encode1based)
        pos0 = pos1 - 1
        alt_norm = alt.upper()
        ref_work = ref

        # Strip leading common characters
        while len(ref_work) > 0 and len(alt_norm) > 0 and ref_work[0] == alt_norm[0]:
            pos0 += 1
            ref_work = ref_work[1:]
            alt_norm = alt_norm[1:]

        del_len = len(ref_work)

        # encodeSpdi
        result2 = encoder.encodeSpdi(chrom, pos0, del_len, alt_norm)

        print(f"  SPDI: {chrom}:{pos0} del={del_len} ins='{alt_norm}'")
        print(f"  encode1based: {result1}")
        print(f"  encodeSpdi:   {result2}")
        print(f"  Match: {result1 == result2}")
        print()

    print("3. Direct SPDI usage examples")
    print("-" * 32)

    spdi_examples = [
        # (description, chr, pos_0based, del_length, ins)
        ("Pure deletion", "1", 1000, 5, ""),
        ("Pure insertion", "1", 1001, 0, "ATCG"),
        ("Substitution", "1", 1002, 3, "GGG"),
        ("Complex indel", "1", 1003, 2, "TTTTTT"),
    ]

    for desc, chrom, pos0, del_len, ins in spdi_examples:
        result = encoder.encodeSpdi(chrom, pos0, del_len, ins)
        print(f"{desc}: {chrom}:{pos0} del={del_len} ins='{ins}' -> {result}")

    print()
    print("4. Advantages of encodeSpdi")
    print("-" * 29)
    print("- Direct SPDI format input (no conversion needed)")
    print("- Clearer semantics for deletions and insertions")
    print("- No need to construct reference sequences")
    print("- Better for programmatic variant generation")
    print("- Follows DRY principle (encode1based now uses encodeSpdi internally)")


if __name__ == "__main__":
    main()
