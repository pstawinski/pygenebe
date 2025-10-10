# entrypoint.py
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
        logging.error(
            f"Not supported anymore, use https://github.com/pstawinski/genebe-cli or let us know if you need this feature."
        )
    elif args.command == "account":
        account_command(args)
    elif args.command == "version":
        version_command(args)
    else:
        logging.error(f"Not recognized command: {args.command}")


if __name__ == "__main__":
    main()
