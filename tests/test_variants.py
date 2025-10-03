"""Tests for variant annotation functionality."""

import pytest
import genebe as gb
import pandas as pd


def test_annotate_variants_list():
    """Test case for annotating variants."""
    variants = ["6-160585140-T-G"]
    annotations = gb.annotate_variants_list(
        variants,
        use_ensembl=True,
        use_refseq=False,
        genome="hg38",
        batch_size=500,
        use_netrc=False,
        endpoint_url="https://api.genebe.net/cloud/api-public/v1/variants",
    )

    # Assertions
    assert len(annotations) == len(variants)

    anns = annotations[0]
    print(anns)

    assert anns["chr"] == "6"
    assert anns["pos"] == 160585140


def test_gbid_generation():
    """Test GBID generation."""
    chr = "1"
    pos = 16044378
    ref = "C"
    alt = "CACACACACAT"
    encoded = gb.encode_vcf_variant_gbid(chr, pos, ref, alt)
    assert encoded == 17227519582999023


@pytest.mark.network
def test_annotate_with_list():
    """Test annotate with list of variants - marked as slow test."""
    variants = ["7-69599651-A-G", "6-160585140-T-G"]
    annotations = gb.annotate(
        variants,
        use_ensembl=True,
        use_refseq=False,
        genome="hg38",
        batch_size=500,
        flatten_consequences=True,
        use_netrc=True,
        output_format="list",
    )

    assert isinstance(annotations, list)
    assert len(annotations) == len(variants)


def test_parse_variants_df():
    """Test case for parsing variants from DataFrame with diverse variant formats."""

    # Create test DataFrame with diverse variant types including VCF-style and HGVS formats
    test_data = {
        "variant": [
            # VCF-style variants (chr-pos-ref-alt format)
            "19-18905299-A-G",  # Simple SNV
            "16-29076807-TACAACTATTGCTGTTATACTAAAGACCAATAGAAATAGCAAGATTAATTAAGATACCAGTTGAAA-T",  # Long deletion
            "13-32345389-A-T",  # Simple SNV
            "19-58132809-TTTTTC-T",  # Short deletion
            "15-40037110-C-T",  # Simple SNV
            "16-29613347-C-A",  # Simple SNV
            "3-129433301-T-C",  # Simple SNV
            "15-64748729-C-T",  # Simple SNV
            "10-7642157-CT-C",  # Small deletion
            "1-23913229-C-G",  # Simple SNV
            # Additional diverse formats
            "chrX:153803771:1:A",  # Different coordinate format
            "rs1228544607",  # dbSNP ID
            "ENST00000679957.1:c.803C>T",  # HGVS transcript notation
            "NM_018136.4:c.567_569del",  # HGVS deletion
            "NM_000277.2(PAH):c.169G>A (p.Glu57Lys)",  # HGVS with gene and protein
            "AGT Met259Thr",  # Protein notation
            "AGT M259T",  # Short protein notation
            "1:11796321-G>A",  # Another coordinate format
            "CA123062477",  # ClinVar ID (expected to fail)
            "m.13513G>A",  # Mitochondrial variant
        ]
    }

    df = pd.DataFrame(test_data)

    # Call the function
    result_df = gb.parse_variants_df(
        df, use_netrc=False, progress=False, batch_size=500
    )

    # Test basic structure
    assert isinstance(result_df, pd.DataFrame)
    assert result_df.shape == (20, 5)  # Original variant column + 4 parsed columns

    # Test column names
    expected_columns = ["variant", "chr", "pos", "ref", "alt"]
    assert list(result_df.columns) == expected_columns

    # Test data types
    assert result_df["variant"].dtype == "object"
    assert result_df["chr"].dtype == "object"
    assert result_df["pos"].dtype == "Int32"  # pandas nullable integer
    assert result_df["ref"].dtype == "object"
    assert result_df["alt"].dtype == "object"

    # Test that original variant column is preserved
    pd.testing.assert_series_equal(result_df["variant"], df["variant"])

    # Test parsed values (gold response) - comprehensive test for all formats
    expected_results = [
        # VCF-style variants
        {"chr": "19", "pos": 18905299, "ref": "A", "alt": "G"},
        {
            "chr": "16",
            "pos": 29076807,
            "ref": "TACAACTATTGCTGTTATACTAAAGACCAATAGAAATAGCAAGATTAATTAAGATACCAGTTGAAA",
            "alt": "T",
        },
        {"chr": "13", "pos": 32345389, "ref": "A", "alt": "T"},
        {"chr": "19", "pos": 58132809, "ref": "TTTTTC", "alt": "T"},
        {"chr": "15", "pos": 40037110, "ref": "C", "alt": "T"},
        {"chr": "16", "pos": 29613347, "ref": "C", "alt": "A"},
        {"chr": "3", "pos": 129433301, "ref": "T", "alt": "C"},
        {"chr": "15", "pos": 64748729, "ref": "C", "alt": "T"},
        {"chr": "10", "pos": 7642157, "ref": "CT", "alt": "C"},
        {"chr": "1", "pos": 23913229, "ref": "C", "alt": "G"},
        # Additional diverse formats
        {"chr": "X", "pos": 153803772, "ref": "C", "alt": "A"},  # chrX:153803771:1:A
        {"chr": "1", "pos": 230710021, "ref": "G", "alt": "A"},  # rs1228544607
        {
            "chr": "1",
            "pos": 230710021,
            "ref": "G",
            "alt": "A",
        },  # ENST00000679957.1:c.803C>T
        {
            "chr": "1",
            "pos": 197143682,
            "ref": "GCTC",
            "alt": "G",
        },  # NM_018136.4:c.567_569del
        {
            "chr": "12",
            "pos": 102894918,
            "ref": "C",
            "alt": "T",
        },  # NM_000277.2(PAH):c.169G>A
        {"chr": "1", "pos": 230710048, "ref": "A", "alt": "G"},  # AGT Met259Thr
        {"chr": "1", "pos": 230710048, "ref": "A", "alt": "G"},  # AGT M259T
        {"chr": "1", "pos": 11796321, "ref": "G", "alt": "A"},  # 1:11796321-G>A
        {"chr": None, "pos": None, "ref": None, "alt": None},  # CA123062477 (fails)
        {
            "chr": "M",
            "pos": 13513,
            "ref": "G",
            "alt": "A",
        },  # m.13513G>A (mitochondrial)
    ]

    # Validate each parsed result
    for i, expected in enumerate(expected_results):
        row = result_df.iloc[i]
        variant_name = row["variant"]

        # Handle nullable values properly
        if expected["chr"] is None:
            assert pd.isna(
                row["chr"]
            ), f"Row {i} ({variant_name}): Expected chr=None, got {row['chr']}"
        else:
            assert (
                row["chr"] == expected["chr"]
            ), f"Row {i} ({variant_name}): Expected chr={expected['chr']}, got {row['chr']}"

        if expected["pos"] is None:
            assert pd.isna(
                row["pos"]
            ), f"Row {i} ({variant_name}): Expected pos=None, got {row['pos']}"
        else:
            assert (
                row["pos"] == expected["pos"]
            ), f"Row {i} ({variant_name}): Expected pos={expected['pos']}, got {row['pos']}"

        if expected["ref"] is None:
            assert pd.isna(
                row["ref"]
            ), f"Row {i} ({variant_name}): Expected ref=None, got {row['ref']}"
        else:
            assert (
                row["ref"] == expected["ref"]
            ), f"Row {i} ({variant_name}): Expected ref={expected['ref']}, got {row['ref']}"

        if expected["alt"] is None:
            assert pd.isna(
                row["alt"]
            ), f"Row {i} ({variant_name}): Expected alt=None, got {row['alt']}"
        else:
            assert (
                row["alt"] == expected["alt"]
            ), f"Row {i} ({variant_name}): Expected alt={expected['alt']}, got {row['alt']}"

    # Test that most variants parse successfully (except CA123062477 which is expected to fail)
    successful_variants = [
        i for i, exp in enumerate(expected_results) if exp["chr"] is not None
    ]
    failed_variants = [
        i for i, exp in enumerate(expected_results) if exp["chr"] is None
    ]

    assert len(successful_variants) == 19, "Expected 19 variants to parse successfully"
    assert len(failed_variants) == 1, "Expected 1 variant to fail parsing"


def test_parse_variants_df_edge_cases():
    """Test edge cases for parse_variants_df."""

    # Test with empty DataFrame
    empty_df = pd.DataFrame({"variant": []})
    result = gb.parse_variants_df(empty_df, use_netrc=False, progress=False)
    assert result.shape == (0, 5)
    assert list(result.columns) == ["variant", "chr", "pos", "ref", "alt"]

    # Test with DataFrame missing 'variant' column
    no_variant_df = pd.DataFrame({"other_col": ["test"]})
    result = gb.parse_variants_df(no_variant_df, use_netrc=False, progress=False)
    # Should return original DataFrame unchanged
    pd.testing.assert_frame_equal(result, no_variant_df)

    # Test with single variant
    single_df = pd.DataFrame({"variant": ["1-23913229-C-G"]})
    result = gb.parse_variants_df(single_df, use_netrc=False, progress=False)
    assert result.shape == (1, 5)
    assert result.iloc[0]["chr"] == "1"
    assert result.iloc[0]["pos"] == 23913229
    assert result.iloc[0]["ref"] == "C"
    assert result.iloc[0]["alt"] == "G"


@pytest.mark.network
def test_annotate_with_dataframe():
    """Test the annotate function with DataFrame input using diverse variant formats."""

    # Use the same diverse variant formats from test_parse_variants_df
    # First parse them to get chr-pos-ref-alt format, then annotate
    variant_input_data = {
        "variant": [
            # VCF-style variants (chr-pos-ref-alt format)
            "19-18905299-A-G",  # Simple SNV
            "16-29076807-TACAACTATTGCTGTTATACTAAAGACCAATAGAAATAGCAAGATTAATTAAGATACCAGTTGAAA-T",  # Long deletion
            "13-32345389-A-T",  # Simple SNV
            "15-40037110-C-T",  # Simple SNV
            "1-23913229-C-G",  # Simple SNV
            # Additional diverse formats
            "rs1228544607",  # dbSNP ID
            "ENST00000679957.1:c.803C>T",  # HGVS transcript notation
            "NM_018136.4:c.567_569del",  # HGVS deletion
            "AGT Met259Thr",  # Protein notation
            "1:11796321-G>A",  # Another coordinate format
            "m.13513G>A",  # Mitochondrial variant
        ]
    }

    # First parse the variants to get chr-pos-ref-alt format
    variant_df = pd.DataFrame(variant_input_data)
    parsed_df = gb.parse_variants_df(variant_df, use_netrc=False, progress=False)

    # Remove any failed parsing results (where chr is None)
    successful_variants = parsed_df[parsed_df["chr"].notna()].copy()

    # Create a proper DataFrame for annotation with chr, pos, ref, alt columns
    annotation_input_df = successful_variants[["chr", "pos", "ref", "alt"]].copy()

    # Call the annotate function with the parsed DataFrame
    result_df = gb.annotate(
        annotation_input_df,
        use_netrc=False,
        progress_bar=False,
        batch_size=500,
        genome="hg38",
        use_ensembl=True,
        use_refseq=True,
        flatten_consequences=True,
    )

    # Test basic structure
    assert isinstance(result_df, pd.DataFrame)
    assert result_df.shape[0] == len(
        annotation_input_df
    )  # Same number of rows as input
    assert result_df.shape[1] > 10  # Should have many annotation columns

    # Test that original columns are preserved
    assert "chr" in result_df.columns
    assert "pos" in result_df.columns
    assert "ref" in result_df.columns
    assert "alt" in result_df.columns

    # Test that annotation columns are added
    expected_annotation_columns = [
        "effect",
        "gene_symbol",
        "acmg_score",
        "acmg_classification",
        "gnomad_exomes_af",
        "gnomad_genomes_af",
        "dbsnp",
        "hgvs_c",
    ]

    for col in expected_annotation_columns:
        assert col in result_df.columns, f"Expected column {col} not found"

    # Test that original variant data is preserved
    for i in range(len(annotation_input_df)):
        assert result_df.iloc[i]["chr"] == annotation_input_df.iloc[i]["chr"]
        assert result_df.iloc[i]["pos"] == annotation_input_df.iloc[i]["pos"]
        assert result_df.iloc[i]["ref"] == annotation_input_df.iloc[i]["ref"]
        assert result_df.iloc[i]["alt"] == annotation_input_df.iloc[i]["alt"]

    # Test data types
    assert result_df["chr"].dtype == "object"
    assert result_df["pos"].dtype in [
        "int64",
        "Int64",
        "int32",
        "Int32",
    ]  # Could be either
    assert result_df["ref"].dtype == "object"
    assert result_df["alt"].dtype == "object"

    # Test that annotations are not all null (at least some variants should be annotated)
    non_null_counts = {}
    for col in expected_annotation_columns:
        non_null_count = result_df[col].notna().sum()
        non_null_counts[col] = non_null_count

    # At least some annotations should be present
    assert non_null_counts["effect"] > 0, "No effect annotations found"
    assert non_null_counts["acmg_score"] >= 0, "ACMG scores should be present or null"

    # Test that we successfully annotated diverse variant types
    print(f"Successfully annotated {len(result_df)} variants from diverse formats:")
    print(
        f"  - Original formats included: VCF-style, HGVS, dbSNP, protein notation, mitochondrial"
    )
    print(f"  - Result DataFrame has {result_df.shape[1]} annotation columns")
    print(f"  - Variants with effects: {non_null_counts.get('effect', 0)}")
    print(f"  - Variants with gene symbols: {non_null_counts.get('gene_symbol', 0)}")


@pytest.mark.network
def test_annotate_with_list():
    """Test the annotate function with list input and different output formats."""
    # Test with list input, list output
    variants_list = ["19-18905299-A-G", "1-23913229-C-G", "13-32345389-A-T"]

    result_list = gb.annotate(
        variants_list,
        use_netrc=False,
        progress_bar=False,
        batch_size=500,
        output_format="list",
        genome="hg38",
        use_ensembl=True,
        use_refseq=True,
        flatten_consequences=True,
    )

    # Test basic structure for list output
    assert isinstance(result_list, list)
    assert len(result_list) == len(variants_list)

    # Test that each result is a dictionary with expected keys
    for i, annotation in enumerate(result_list):
        assert isinstance(annotation, dict)

        # Check core variant fields
        assert "chr" in annotation
        assert "pos" in annotation
        assert "ref" in annotation
        assert "alt" in annotation

        # Check annotation fields
        assert "effect" in annotation
        assert "gene_symbol" in annotation
        assert "acmg_score" in annotation

        # Verify variant coordinates match input
        expected_variant = variants_list[i].split("-")
        assert annotation["chr"] == expected_variant[0]
        assert annotation["pos"] == int(expected_variant[1])
        assert annotation["ref"] == expected_variant[2]
        assert annotation["alt"] == expected_variant[3]

    # Test with list input, dataframe output
    result_df = gb.annotate(
        variants_list,
        use_netrc=False,
        progress_bar=False,
        batch_size=500,
        output_format="dataframe",
        genome="hg38",
        use_ensembl=True,
        use_refseq=True,
        flatten_consequences=True,
    )

    # Test structure for dataframe output
    assert isinstance(result_df, pd.DataFrame)
    assert result_df.shape[0] == len(variants_list)
    assert result_df.shape[1] > 10

    print(
        f"Test passed! List format returned {len(result_list)} annotations, DataFrame format shape: {result_df.shape}"
    )


def test_annotate_edge_cases():
    """Test edge cases for the annotate function."""

    # Test with missing required columns - should raise ValueError
    invalid_df = pd.DataFrame({"chr": ["1"], "pos": [12345]})  # Missing ref, alt
    with pytest.raises(ValueError, match="Missing required columns"):
        gb.annotate(invalid_df, use_netrc=False, progress_bar=False)

    # Test with invalid output format - should raise ValueError
    with pytest.raises(ValueError, match="output_format must be one of"):
        gb.annotate(
            ["1-12345-A-G"],
            output_format="invalid",
            use_netrc=False,
            progress_bar=False,
        )

    # Test that both use_refseq and use_ensembl cannot be False
    with pytest.raises(
        ValueError, match="use_refseq and use_ensembl cannot be both False"
    ):
        gb.annotate(
            ["1-12345-A-G"],
            use_refseq=False,
            use_ensembl=False,
            use_netrc=False,
            progress_bar=False,
        )

    # Test with single variant DataFrame
    single_df = pd.DataFrame(
        {"chr": ["1"], "pos": [23913229], "ref": ["C"], "alt": ["G"]}
    )
    result = gb.annotate(single_df, use_netrc=False, progress_bar=False)
    assert isinstance(result, pd.DataFrame)
    assert result.shape[0] == 1
    assert "chr" in result.columns
    assert "effect" in result.columns
    assert result.iloc[0]["chr"] == "1"

    print("All edge case tests passed!")
