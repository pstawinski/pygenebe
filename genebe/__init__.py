# may be empty

from .gbid import (
    encode_vcf_position_gbid,
    encode_vcf_variant_gbid,
)

from .transcriptid import (
    encode_transcript_id,
    decode_transcript_id,
)

from .client import (
    annotate_variants_list,
    annotate_variants_list_to_dataframe,
    parse_variants,
    parse_variants_df,
    annotate,
    annotate_dataframe_variants,
    whoami,
    lift_over_variants,
    lift_over_variants_df,
)
from .vcf_simple_annotator import annotate_vcf
from .json_simple_annotator import annotate_json
from .version import __version__
