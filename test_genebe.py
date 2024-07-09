from genebe import annotate

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


if __name__ == "__main__":
    test_annotate_with_dataframe()


if __name__ == "__main__":
    test_annotate_with_list()
    test_annotate_with_dataframe()
    print("All tests passed!")
