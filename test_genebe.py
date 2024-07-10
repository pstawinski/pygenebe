from genebe import annotate, parse_hgvs, parse_spdi

import pandas as pd


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


def test_parse_hgvs():
    # Test parse_hgvs with a list of HGVS variants
    hgvs_variants = ["NM_000277.2:c.1A>G", "NM_000277.2:c.2T>C"]
    annotations = parse_hgvs(
        hgvs_variants,
        batch_size=500,
        username=None,
        api_key=None,
        use_netrc=True,
        endpoint_url="https://api.genebe.net/cloud/api-public/v1/hgvs",
    )
    print("HGVS Annotations:")
    print(annotations)
    assert isinstance(annotations, list)
    assert len(annotations) == len(
        hgvs_variants
    )  # Ensure number of annotations matches number of variants
    for annotation in annotations:
        assert isinstance(annotation, str)  # Ensure each annotation is a string


def test_parse_spdi():
    # Test parse_spdi with a list of SPDI variants
    spdi_variants = ["chrX:153803771:1:A"]
    annotations = parse_spdi(
        spdi_variants,
        genome="hg38",
        batch_size=500,
        username=None,
        api_key=None,
        use_netrc=True,
        endpoint_url="https://api.genebe.net/cloud/api-public/v1/spdi",
    )
    print("SPDI Annotations:")
    print(annotations)
    assert isinstance(annotations, list)
    assert len(annotations) == len(
        spdi_variants
    )  # Ensure number of annotations matches number of variants
    for annotation in annotations:
        assert isinstance(annotation, str)  # Ensure each annotation is a string


if __name__ == "__main__":
    test_annotate_with_list()
    test_annotate_with_list_custom_annotations()
    test_annotate_with_dataframe()
    test_parse_hgvs()
    test_parse_spdi()
    print("All tests passed!")
