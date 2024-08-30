# entrypoint.py

from .vcf_simple_annotator import annotate_vcf
from .json_simple_annotator import annotate_json
from .version import __version__
from .client import whoami
import argparse
import logging


def account_command(args):
    # Add your login logic here
    print(f"Logging in with username: {args.username} and API key: {args.api_key}")
    result = whoami(
        api_key=args.api_key,
        use_netrc=args.use_netrc,
        username=args.username,
        endpoint_url=args.endpoint_url,
    )
    print(result)


def version_command(args):
    # Add your login logic here
    print(f"GeneBe, version {__version__}")


def main():
    parser = argparse.ArgumentParser(
        description="GeneBe client, v" + __version__ + " :: https://genebe.net",
    )

    subparsers = parser.add_subparsers(title="Commands", dest="command")
    subparsers.required = True
    # Subparser for the 'annotate' command
    annotate_parser = subparsers.add_parser("annotate", help="Annotate VCF or JSON")

    annotate_parser.add_argument(
        "--input", required=True, help="Path to the input file"
    )
    annotate_parser.add_argument(
        "--output", required=True, help="Path to the output file"
    )
    annotate_parser.add_argument(
        "--input_type",
        default="vcf",
        help="Input and output file type. Possible values: {vcf, json} (default: vcf)",
    )
    annotate_parser.add_argument(
        "--genome", default="hg38", help="Genome version (default: hg38)"
    )
    annotate_parser.add_argument(
        "--omit_ensembl",
        help="Use Ensembl data for annotation (True/False). Default is True.",
        action="store_true",
        default=False,
    )
    annotate_parser.add_argument(
        "--omit_refseq",
        help="Use RefSeq data for annotation (True/False). Default is True.",
        action="store_true",
        default=False,
    )
    annotate_parser.add_argument(
        "--flatten_consequences",
        action="store_true",
        help="Flatten consequences in the output",
    )
    annotate_parser.add_argument(
        "--batch_size",
        type=int,
        default=500,
        help="Batch size for API requests (default: 500)",
    )
    annotate_parser.add_argument("--username", help="Username for API authentication")
    annotate_parser.add_argument("--api_key", help="API key for authentication")
    annotate_parser.add_argument(
        "--use_netrc", action="store_true", help="Use .netrc file for authentication"
    )
    annotate_parser.add_argument(
        "--endpoint_url",
        default="https://api.genebe.net/cloud/api-public/v1/variants",
        help="API endpoint URL (default: https://api.genebe.net/cloud/api-public/v1/variants)",
    )
    annotate_parser.add_argument(
        "--omit_acmg",
        action="store_true",
        help="Don't add ACMG scores in the output",
        default=False,
    )
    annotate_parser.add_argument(
        "--omit_csq",
        action="store_true",
        help="Don't add consequences in the output",
        default=False,
    )
    annotate_parser.add_argument(
        "--omit_basic",
        action="store_true",
        help="Don't add basic annotations (GnomAD frequencies etc) in the output",
        default=False,
    )
    annotate_parser.add_argument(
        "--omit_advanced",
        action="store_true",
        help="Don't add advanced annotations (ClinVar frequencies etc) in the output",
        default=False,
    )
    annotate_parser.add_argument(
        "--omit_normalization",
        action="store_true",
        help="Don't normalize variants. Use only if you are sure they are normalized already",
        default=False,
    )
    annotate_parser.add_argument(
        "--progress",
        action="store_true",
        help="Show progress bar",
        default=False,
    )

    # Subparser for the 'account' command
    account_parser = subparsers.add_parser(
        "account", help="Account information in GeneBe"
    )
    account_parser.add_argument("--username", help="Username")
    account_parser.add_argument("--api_key", help="API key")
    account_parser.add_argument(
        "--use_netrc", action="store_true", help="Use .netrc file for authentication"
    )
    account_parser.add_argument(
        "--endpoint_url",
        default="https://api.genebe.net/cloud/api-public/v1/whoami",
        help="API endpoint URL (default: https://api.genebe.net/cloud/api-public/v1/whoami)",
    )

    # Subparser for the 'version' command
    account_parser = subparsers.add_parser("version", help="Display client version")

    args = parser.parse_args()

    # Call the appropriate function based on the selected command
    if args.command == "annotate":
        # Call the annotate_vcf function with the provided arguments
        type = args.input_type
        if type == "vcf":
            annotate_vcf(
                input_vcf_path=args.input,
                output_vcf_path=args.output,
                genome=args.genome,
                use_ensembl=not args.omit_refseq,
                use_refseq=not args.omit_refseq,
                flatten_consequences=args.flatten_consequences,
                batch_size=args.batch_size,
                username=args.username,
                api_key=args.api_key,
                use_netrc=args.use_netrc,
                endpoint_url=args.endpoint_url,
                omit_acmg=args.omit_acmg,
                omit_csq=args.omit_csq,
                omit_basic=args.omit_basic,
                omit_advanced=args.omit_advanced,
                omit_normalization=args.omit_normalization,
                progress_bar=args.progress,
            )
        elif type == "json":
            annotate_json(
                input_json_path=args.input,
                output_json_path=args.output,
                genome=args.genome,
                use_ensembl=not args.omit_ensembl,
                use_refseq=not args.omit_refseq,
                flatten_consequences=args.flatten_consequences,
                batch_size=args.batch_size,
                username=args.username,
                api_key=args.api_key,
                use_netrc=args.use_netrc,
                endpoint_url=args.endpoint_url,
                omit_acmg=args.omit_acmg,
                omit_csq=args.omit_csq,
                omit_basic=args.omit_basic,
                omit_advanced=args.omit_advanced,
                omit_normalization=args.omit_normalization,
                progress_bar=args.progress,
            )
        else:
            logging.error(
                f"Invalid input type given: {type}. Check help for supported file types."
            )
    elif args.command == "account":
        account_command(args)
    elif args.command == "version":
        version_command(args)
    else:
        logging.error(f"Not recognized command: {args.command}")


if __name__ == "__main__":
    main()
