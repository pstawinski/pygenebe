from genebe import (
    annotate,
    parse_variants,
    lift_over_variants,
    annotate_dataframe_variants,
    encode_vcf_variant_gbid,
)

import pandas as pd


def test_gbid_generation():
    chr = "1"
    pos = 16044378
    ref = "C"
    alt = "CACACACACAT"
    encoded = encode_vcf_variant_gbid(chr, pos, ref, alt)
    assert encoded == 17227519582999023


def test_annotate_with_list():
    # Test annotate with list of variants
    variants = ["7-69599651-A-G", "6-160585140-T-G"]
    annotations = annotate(
        variants,
        use_ensembl=True,
        use_refseq=False,
        genome="hg38",
        batch_size=500,
        flatten_consequences=True,
        use_netrc=True,
        output_format="list",
    )
    print("Annotations (list format):")
    print(annotations)
    assert isinstance(annotations, list)
    assert len(annotations) == len(
        variants
    )  # Ensure number of annotations matches number of variants


def test_annotate_with_list_custom_annotations():
    # Test annotate with list of variants
    variants = ["7-69599651-A-G", "6-160585140-T-G"]
    annotations = annotate(
        variants,
        use_ensembl=True,
        use_refseq=False,
        genome="hg38",
        batch_size=500,
        flatten_consequences=True,
        use_netrc=True,
        output_format="list",
        custom_annotations=["AC_gnomad4", "AN_gnomad4"],
    )
    print("Annotations (list format) with custom annotations:")
    print(annotations)
    assert isinstance(annotations, list)
    assert len(annotations) == len(variants)


def test_annotate_with_dataframe():

    # Test annotate with pandas DataFrame
    data = {
        "chr": ["7", "6"],
        "pos": [69599651, 160585140],
        "ref": ["A", "T"],
        "alt": ["G", "G"],
    }
    df = pd.DataFrame(data)

    annotations_df = annotate(
        df,
        use_ensembl=True,
        use_refseq=False,
        genome="hg38",
        batch_size=500,
        flatten_consequences=True,
        use_netrc=True,
        output_format="dataframe",
    )
    print("Annotations (DataFrame format):")
    print(annotations_df)

    assert isinstance(annotations_df, pd.DataFrame)
    assert set(["chr", "pos", "ref", "alt"]).issubset(set(annotations_df.columns))
    assert (
        annotations_df.shape[0] == df.shape[0]
    )  # Ensure number of rows matches number of variants


def test_lift_over():
    variants = ["chr6-160585140-T-G", "chrX-200-C-G", "invalid_input_missing_dash"]

    vars = lift_over_variants(
        variants,
        "hg19",
        "hg38",
        endpoint_url="http://localhost:7180/cloud/api-public/v1/liftover",
    )

    print("Lifted over variants:")
    print(vars)
    print(len(vars))
    print(len(variants))
    assert len(vars) == len(variants)


def test_parse_variants():
    # Test parse_spdi with a list of SPDI variants
    variants = [
        "chrX:153803771:1:A",
        "22 28695868 AG A",
        "22-28695869--G",
        "22-28695869-G-",
        "NM_000277.2:c.1A>G",
        "NM_000277.2:c.2T>C",
        "AGT M259T",
        "rs1228544607",
        "wadliwy",
    ]
    annotations = parse_variants(
        variants,
        genome="hg38",
        batch_size=500,
        username=None,
        api_key=None,
        use_netrc=True,
        endpoint_url="https://api.genebe.net/cloud/api-public/v1/convert",
    )
    print("Variants:")
    print(annotations)
    assert isinstance(annotations, list)
    assert len(annotations) == len(
        variants
    )  # Ensure number of annotations matches number of variants
    for annotation in annotations:
        assert isinstance(annotation, str)  # Ensure each annotation is a string


def test_annotate_variants_list_hg19():
    # Test case
    variants = ["1-12979845-T-C"]
    annotations = annotate(
        variants,
        use_ensembl=True,
        use_refseq=False,
        genome="hg19",
        batch_size=500,
        use_netrc=False,
        endpoint_url="https://api.genebe.net/cloud/api-public/v1/variants",
    )

    # Assertions
    assert len(annotations) == len(variants)

    anns = annotations[0]
    print(anns)

    assert anns["chr"] == "1"
    assert anns["pos"] == 12920024


def test_annotate_variants_list_hg19_impossible_liftover():
    # Test case
    variants = ["1-12979845-T-C", "1-13036247-T-A"]
    annotations = annotate(
        variants,
        use_ensembl=True,
        use_refseq=False,
        genome="hg19",
        batch_size=500,
        use_netrc=False,
        endpoint_url="https://api.genebe.net/cloud/api-public/v1/variants",
        output_format="dataframe",
    )

    # Assertions
    assert len(annotations) == len(variants)

    # anns = annotations[0]
    print(annotations)


def test_annotate_variants_list_hg19_impossible_liftover_old_inteface():
    # Test case
    variants = ["1-12979845-T-C", "1-13036247-T-A"]

    # Split each variant into its components and create a DataFrame
    data = [variant.split("-") for variant in variants]
    df = pd.DataFrame(data, columns=["chr", "pos", "ref", "alt"])

    annotations = annotate_dataframe_variants(
        df,
        use_ensembl=True,
        use_refseq=False,
        genome="hg19",
        batch_size=500,
        use_netrc=False,
        endpoint_url="https://api.genebe.net/cloud/api-public/v1/variants",
    )

    # Assertions
    assert len(annotations) == len(variants)

    # anns = annotations[0]
    print(annotations)


if __name__ == "__main__":
    test_gbid_generation()
    test_annotate_with_list()
    test_annotate_with_list_custom_annotations()
    test_annotate_with_dataframe()
    test_parse_variants()
    test_lift_over()
    test_annotate_variants_list_hg19_impossible_liftover()
    print("All tests passed!")
