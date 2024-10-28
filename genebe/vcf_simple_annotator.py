import logging
from itertools import islice
from .client import annotate
from typing import Optional
from tqdm import tqdm

# Configure logging to print INFO level and above
logging.basicConfig(level=logging.INFO)

try:
    # try with a fast c-implementation ...
    import cyvcf2 as cyvcf2
except ImportError:
    # ... otherwise fallback python implementation
    logging.error(
        'No cyvcf2 module installed. Some features are not available. Install it executing "pip install cyvcf2".'
    )


def annotate_vcf(
    input_vcf_path: str,
    output_vcf_path: str,
    genome: str = "hg38",
    use_ensembl: bool = True,
    use_refseq: bool = True,
    flatten_consequences: bool = True,
    batch_size: int = 500,
    username: Optional[str] = None,
    api_key: Optional[str] = None,
    use_netrc: bool = True,
    endpoint_url: str = "https://api.genebe.net/cloud/api-public/v1/variants",
    omit_acmg: bool = False,
    omit_csq: bool = False,
    omit_basic: bool = False,
    omit_advanced: bool = False,
    omit_normalization: bool = False,
    progress_bar: bool = False,
):
    """
    Annotates a VCF file using the annotate_variants_list function.

    Args:
        input_vcf_path (str): Path to the input VCF file.
        output_vcf_path (str): Path to the output annotated VCF file.
        genome (str, optional): The genome version for annotation (e.g., 'hg38').
            Defaults to 'hg38'.
        use_ensembl (bool, optional): Whether to use Ensembl for annotation.
            Defaults to True.
        use_refseq (bool, optional): Whether to use RefSeq for annotation.
            Defaults to True.
        flatten_consequences (bool, optional): If set to False, return consequences as a list of objects.
            If set to True, only the most important consequence is returned in a flat form.
            Defaults to True.
        batch_size (int, optional): The size of each batch for processing variants.
            Defaults to 100. Must be smaller or equal 1000.
        username (str, optional): The username for authentication.
            Defaults to None.
        api_key (str, optional): The API key for authentication.
            Defaults to None.
        use_netrc (bool, optional): Whether to use credentials from the user's
            .netrc file for authentication. Defaults to True.
        endpoint_url (str, optional): The API endpoint for variant annotation.
            Defaults to 'https://api.genebe.net/cloud/api-public/v1/variants'.
        omit_acmg (bool, optional): Don't add ACMG scores in the output. Defaults to False.
        omit_csq (bool, optional): Don't add consequences in the output. Defaults to False.
        omit_basic (bool, optional): Don't add basic annotations (GnomAD frequencies etc) in the output. Defaults to False.
        omit_advanced (bool, optional): Don't add advanced annotations (ClinVar etc) in the output. Defaults to False.
        omit_normalization (bool, optional): Don't normalize variants. Use only if you are sure they are normalized already. Defaults to False.
        progress (bool, optional): Show progress bar
    """

    logging.info(
        f"Welcome. Annotating VCF file for genome {genome}. Ensure that VCF file has biallelic form."
    )

    # Open VCF file
    vcf_reader = cyvcf2.VCF(input_vcf_path)

    # Add header to the VCF file
    vcf_reader.add_info_to_header(
        {
            "ID": "acmg_score",
            "Number": "1",
            "Type": "Integer",
            "Description": "ACMG Score by GeneBe",
        }
    )
    vcf_reader.add_info_to_header(
        {
            "ID": "acmg_criteria",
            "Number": ".",
            "Type": "String",
            "Description": "ACMG Score by GeneBe",
        }
    )
    vcf_reader.add_info_to_header(
        {
            "ID": "gnomad_exomes_af",
            "Number": "1",
            "Type": "Float",
            "Description": "GnomAD exomes Allele Frequency, by GeneBe",
        }
    )
    vcf_reader.add_info_to_header(
        {
            "ID": "gnomad_genomes_af",
            "Number": "1",
            "Type": "Float",
            "Description": "GnomAD genomes Allele Frequency, by GeneBe",
        }
    )

    # Open the output VCF file for writing
    vcf_writer = cyvcf2.Writer(output_vcf_path, vcf_reader)

    # Batch size for reading variants, if greater than 1000 limit to 1000
    if batch_size > 1000:
        batch_size = 1000
        logging.warning("Batch size is limited to 1000")

    if batch_size < 10:
        batch_size = 10
        logging.warning("Batch size must be greater than 10")

    with tqdm() as pbar:
        # Iterate over variants in batches
        while True:
            # Read a batch of variants
            batch = list(islice(vcf_reader, batch_size))
            if not batch:
                break  # Break if there are no more variants

            # identify variant.ALT that has more than 1 element
            for variant in batch:
                if len(variant.ALT) > 1:
                    raise ValueError(
                        f"Variant {variant.CHROM}-{variant.POS}-{variant.REF}-{variant.ALT} has more than 1 ALT. You have to split multi-allelic variants before annotating for example using bcftools: `bcftools norm -m -any input.vcf -o output_split.vcf`."
                    )

            # Process the batch
            variants_batch = [
                f"{variant.CHROM}-{variant.POS}-{variant.REF}-{variant.ALT[0]}"
                for variant in batch
            ]
            annotated_variants = annotate(
                variants_batch,
                genome=genome,
                use_ensembl=use_ensembl,
                use_refseq=use_refseq,
                flatten_consequences=flatten_consequences,
                batch_size=batch_size,
                username=username,
                api_key=api_key,
                use_netrc=use_netrc,
                endpoint_url=endpoint_url,
                progress_bar=False,
                omit_acmg=omit_acmg,
                omit_csq=omit_csq,
                omit_basic=omit_basic,
                omit_advanced=omit_advanced,
                omit_normalization=omit_normalization,
            )

            # Add ACMG scores to the variants and write to the output VCF file
            for variant, annotated_variant in zip(batch, annotated_variants):
                fields_to_copy = [
                    "acmg_score",
                    "gnomad_exomes_af",
                    "gnomad_genomes_af",
                    "acmg_criteria",
                    #                    "consequences_ensembl",
                    #                    "consequences_refseq",
                ]

                for field in fields_to_copy:
                    value = annotated_variant.get(field)
                    if value is not None:
                        variant.INFO[field] = value

                vcf_writer.write_record(variant)
            pbar.update(len(batch))

    # Close the VCF writer
    vcf_writer.close()
