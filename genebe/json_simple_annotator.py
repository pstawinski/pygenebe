from .client import annotate
from typing import Optional
import jsonlines
from tqdm import tqdm


def annotate_json(
    input_json_path: str,
    output_json_path: str,
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
    Annotates a JSON file using the annotate_variants_list function.

    Args:
        input_json_path (str): Path to the input json lines file.
        output_json_path (str): Path to the output annotated json lines file.
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

    # Batch size for reading variants
    batch_size = 1000

    # Open the output VCF file for writing
    with jsonlines.open(input_json_path, "r") as infile, jsonlines.open(
        output_json_path, "w"
    ) as outfile:
        with tqdm() as pbar:
            batch = []
            for json_object in infile:
                # Add the object to the current batch
                batch.append(json_object)

                # If the batch size is reached, process and write to the output file
                if len(batch) == batch_size:
                    variants_batch = [
                        f"{variant.get('chr')}-{variant.get('pos')}-{variant.get('ref')}-{variant.get('alt')}"
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

                    for input, annotated_variant in zip(batch, annotated_variants):
                        fields_to_copy = [
                            "acmg_score",
                            "gnomad_exomes_af",
                            "gnomad_genomes_af",
                            "consequences_ensembl",
                            "consequences_refseq",
                        ]

                        for field in fields_to_copy:
                            value = annotated_variant.get(field)
                            if value is not None:
                                input[field] = value

                    outfile.write_all(batch)
                    pbar.update(len(batch))
                    # Clear the batch for the next set of data
                    batch = []

        # Process any remaining items in the last batch
        if batch:
            variants_batch = [
                f"{variant.get('chr')}-{variant.get('pos')}-{variant.get('ref')}-{variant.get('alt')}"
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

            for input, annotated_variant in zip(batch, annotated_variants):
                fields_to_copy = [
                    "acmg_score",
                    "gnomad_exomes_af",
                    "gnomad_genomes_af",
                    "consequences_ensembl",
                    "consequences_refseq",
                ]

                for field in fields_to_copy:
                    value = annotated_variant.get(field)
                    if value is not None:
                        input[field] = value

            outfile.write_all(batch)
