from typing import List
from genebe import (
    annotate,
    parse_variants,
    lift_over_variants,
    lift_over_variants_df,
    annotate_dataframe_variants,
    encode_vcf_variant_gbid,
    encode_vcf_position_gbid,
)

import pandas as pd
import json

# This is just to test single bugs or custom features; not a full test suite.
# To be run by hand in dev


def test_custom():
    list = annotate(
        ["chr10-27444157-G-C"],
        use_ensembl=True,
        use_refseq=False,
        genome="hg19",
        batch_size=500,
        use_netrc=True,
        flatten_consequences=True,
        output_format="list",
        endpoint_url="http://localhost:7180/cloud/api-public/v1/variants",
        all_genes=True,
    )
    # print list as a json
    print("Annotations (list format):")
    print(json.dumps(list, indent=2))


if __name__ == "__main__":
    test_custom()
    print("All tests passed!")
