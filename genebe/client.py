import pandas as pd
import requests
from requests.auth import HTTPBasicAuth
from typing import Dict, List
import re
import json
import logging
import netrc
from tinynetrc import Netrc
from urllib.parse import urlparse
from importlib import metadata
from .version import __version__
from pathlib import Path


try:
    from tqdm import tqdm
except ImportError:

    def tqdm(iterator, *args, **kwargs):
        return iterator


# Set up the User-Agent header
user_agent = f"GeneBe/{__version__} (requests)"
apikey_username_warning_showed = False


class TooManyRequestsError(Exception):
    pass


class WrongServerResponseError(Exception):
    pass


def _handle_too_many_requests():
    raise TooManyRequestsError(
        "Too many requests (429 error). Please provide a username and API key. "
        "More information at https://genebe.net/api"
    )


def _parse_variant_string(s):
    if len(s) > 0:
        s = s.strip()

    if len(s) == 0:
        return None

    pattern = re.compile(r"^[^-]+-\d+-[ACGT]*-[ACGT]*$")
    if pattern.match(s):
        chr_, pos, ref, alt = s.split("-")
        pos = int(pos)
        return {"chr": chr_, "pos": pos, "ref": ref, "alt": alt}
    else:
        print(f"Error: Invalid format in [{s}]")
        return None


def _save_credentials_to_netrc(machine, login, password):
    try:
        fpath = Path().home().joinpath(".netrc")
        fpath.touch(exist_ok=True)
        netrc = Netrc(str(fpath))

        netrc[machine]["login"] = login
        netrc[machine]["password"] = password
        netrc.save()
        print(
            "Saved credentials to .netrc, you don't need to provide username/api_key anymore."
        )
    except Exception as e:
        print(f"An error occurred while saving credentials to {fpath}: {e}")


def _get_machine_name_from_endpoint(endpoint):
    """
    Extracts the machine name from an endpoint URL.

    Args:
        endpoint (str): The endpoint URL.

    Returns:
        str: The machine name extracted from the endpoint.
    """
    parsed_url = urlparse(endpoint)
    return parsed_url.netloc


def _read_netrc_credentials(endpoint):
    """
    Reads username and password from the .netrc file based on the machine name derived from the endpoint.

    Args:
        endpoint (str): The endpoint URL.

    Returns:
        tuple: A tuple containing (username, account, password). Returns (None, None, None) if no
        matching entry is found in the .netrc file.
    """
    machine_name = _get_machine_name_from_endpoint(endpoint)

    try:
        # Locate the .netrc file in the user's home directory
        netrc_file_path = netrc.netrc()

        # Get the authentication information for the specified machine
        login_info = netrc_file_path.authenticators(machine_name)

        if login_info:
            username, account, password = login_info
            return username, account, password
        else:
            return None, None, None

    except FileNotFoundError:
        print("Error: .netrc file not found.")
        return None, None, None


# Function to annotate variants
def annotate_variants_list(
    variants: List[str],
    genome: str = "hg38",
    use_ensembl: bool = True,
    use_refseq: bool = True,
    flatten_consequences: bool = True,
    batch_size: int = 500,
    username: str = None,
    api_key: str = None,
    use_netrc: bool = True,
    endpoint_url: str = "https://api.genebe.net/cloud/api-public/v1/variants",
    progress_bar: bool = True,
    omit_acmg: bool = False,
    omit_csq: bool = False,
    omit_basic: bool = False,
    omit_advanced: bool = False,
    omit_normalization: bool = False,
    annotator: str = "snpeff",
) -> List[Dict[str, object]]:
    """
    Annotates a list of genetic variants.

    Args:
        variants (List[str]): A list of genetic variants to be annotated.
            Format: chr-pos-ref-alt. Look at examples below.
        use_ensembl (bool, optional): Whether to use Ensembl for annotation.
            Defaults to True.
        use_refseq (bool, optional): Whether to use RefSeq for annotation.
            Defaults to True.
        genome (str, optional): The genome version for annotation (e.g., 'hg38').
            Defaults to 'hg38'.
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
        progress_bar (bool, optional): Show progress bar.
            Defaults to True.
        omit_acmg (bool, optional): Don't add ACMG scores in the output. Defaults to False.
        omit_csq (bool, optional): Don't add consequences in the output. Defaults to False.
        omit_basic (bool, optional): Don't add basic annotations (GnomAD frequencies etc) in the output. Defaults to False.
        omit_advanced (bool, optional): Don't add advanced annotations (ClinVar etc) in the output. Defaults to False.
        omit_normalization (bool, optional): Don't normalize variants. Use only if you are sure they are normalized already. Defaults to False.
        annotator (str, optional): Which VEP implementation to use.
            Defaults to snpeff.

    Returns:
        List[Dict[str, object]]: A list of dictionaries containing annotation information
        for each variant. The dictionary structure may vary. Check the current documentation
        on https://genebe.net/about/api

    Example:
        >>> variants = ["7-69599651-A-G", "6-160585140-T-G"]
        >>> annotations = annotate_variants_list(variants, use_ensembl=True,
        ...                                      use_refseq=False, genome='hg38',
        ...                                      batch_size=500, username="user123@example.com",
        ...                                      api_key="apikey456", use_netrc=False,
        ...                                      endpoint_url='https://api.genebe.net/cloud/api-public/v1/variants')
        >>> print(annotations)
        [{'chr': '7', 'pos':69599651 (...) }]

    Note:
        - The number of the elements in returned list is always equal to the number of queries.
    """

    if (use_refseq != True) and (use_ensembl != True):
        raise ValueError("use_refseq and use_ensembl cannot be both False")

    if (username is not None) and (api_key is not None):
        auth = (username, api_key)
        if use_netrc:
            _save_credentials_to_netrc(endpoint_url, username, api_key)
    elif use_netrc:
        username, _, api_key = _read_netrc_credentials(endpoint_url)
        if (username is not None) and (api_key is not None):
            auth = (username, api_key)
        else:
            auth = None
    else:
        auth = None

    # if username:
    if api_key == None:
        global apikey_username_warning_showed
        if not apikey_username_warning_showed:
            # show it only once
            logging.warning(
                f"You are not logged in to GeneBe. We recommend using username and api_key arguments. Find out more on https://genebe.net/api ."
            )
            apikey_username_warning_showed = True
    else:
        logging.warning(f"you will log in {username}")

    # input data validation
    # Convert the list of strings to list of dictionaries
    dict_list = [_parse_variant_string(s) for s in variants]
    logging.debug("I will query for " + json.dumps(dict_list))

    annotated_data = []

    ranges = range(0, len(dict_list), batch_size)
    iterator = tqdm(ranges) if progress_bar else ranges

    for i in iterator:
        # Prepare data for API request
        api_data = dict_list[i : i + batch_size]

        logging.debug("Querying for " + json.dumps(api_data))

        params = {"genome": genome}
        if use_refseq == False:
            params["useRefseq"] = use_refseq
        if use_ensembl == False:
            params["useEnsembl"] = use_ensembl
        if omit_acmg:
            params["omitAcmg"] = omit_acmg
        if omit_csq:
            params["omitCsq"] = omit_csq
        if omit_basic:
            params["omitBasic"] = omit_basic
        if omit_advanced:
            params["omitAdvanced"] = omit_advanced
        if omit_normalization:
            params["omitNormalization"] = omit_normalization
        if annotator:
            params["annotator"] = annotator

        # Make API request
        response = requests.post(
            endpoint_url,
            params=params,
            json=api_data,
            headers={
                "Accept": "application/json",
                "Content-Type": "application/json",
                "User-Agent": user_agent,
            },
            auth=auth,
        )

        # Check if request was successful
        if response.status_code == 200:
            api_results_raw = response.json()
            api_results = [element for element in api_results_raw["variants"]]

            logging.debug("Backend result is " + json.dumps(api_results))

            # Append API results to annotated_data
            annotated_data.extend(api_results)
        elif response.status_code == 429:
            _handle_too_many_requests()
        else:
            logging.error(
                f"Got response with code {response.status_code} with body "
                + response.text
            )
            raise WrongServerResponseError(
                f"Got response with code {response.status_code} with body "
                + response.text
            )

    if flatten_consequences:
        for item in annotated_data:
            consequences_fields = None
            transcript = item.get("transcript")
            if transcript:  # check the ACMG chosen transcript
                if transcript.startswith("E"):
                    consequences_fields = "consequences_ensembl"
                else:
                    consequences_fields = "consequences_refseq"

                consequences_refseq = item[consequences_fields]
                if consequences_refseq:  # Check if the list is not empty
                    first_consequence = consequences_refseq[0]
                    for key in ["gene_symbol", "gene_hgnc_id", "transcript", "hgvs_c"]:
                        item[key] = first_consequence.get(key, None)
                    effects_list = first_consequence.get("consequences", None)
                    item["consequences"] = (
                        ",".join(
                            str(item) if item is not None else "None"
                            for item in effects_list
                        )
                        if effects_list is not None
                        else ""
                    )
            # Remove the 'consequences_*' fields, if exists
            item.pop("consequences_ensembl", None)
            item.pop("consequences_refseq", None)

    return annotated_data


def annotate_variants_list_to_dataframe(
    variants: List[str],
    use_ensembl: bool = True,
    use_refseq: bool = True,
    genome: str = "hg38",
    batch_size: int = 100,
    flatten_consequences: bool = True,
    username: str = None,
    api_key: str = None,
    use_netrc: bool = True,
    endpoint_url: str = "https://api.genebe.net/cloud/api-public/v1/variants",
) -> pd.DataFrame:
    # Call the existing function
    result_list = annotate_variants_list(
        variants=variants,
        use_ensembl=use_ensembl,
        use_refseq=use_refseq,
        flatten_consequences=flatten_consequences,
        genome=genome,
        batch_size=batch_size,
        username=username,
        api_key=api_key,
        use_netrc=use_netrc,
        endpoint_url=endpoint_url,
    )

    # Convert the list of dictionaries to a Pandas DataFrame
    df = pd.DataFrame(result_list)
    return df


class DnaChange:
    def __init__(self, chr, pos, ref, alt):
        self.chr = chr
        self.pos = int(pos)
        self.ref = ref
        self.alt = alt


def lift_over_variants(
    variants: List[str],
    from_genome: str = "hg19",
    dest_genome: str = "hg38",
    username: str = None,
    api_key: str = None,
    use_netrc: bool = True,
    endpoint_url: str = "https://api.genebe.net/cloud/api-public/v1/liftover",
) -> List[str]:
    """
    Lifts over a list of genetic variants between human genome versions.

    Args:
        variants (List[str]): A list of genetic variants represented as strings (chr-pos-ref-alt).
        from_genome (str): Source human genome version, one of hg19, hg38, t2t.
        dest_genome (str): Destination human genome version, one of hg19, hg38, t2t.
        username (str, optional): The username for authentication.
            Defaults to None.
        api_key (str, optional): The API key for authentication.
            Defaults to None.
        use_netrc (bool, optional): Whether to use credentials from the user's
            .netrc file for authentication. Defaults to True.
        endpoint_url (str, optional): The API endpoint.
            Defaults to 'https://api.genebe.net/cloud/api-public/v1/liftover'.



    Returns:
        List[str]: A list of lifted variants.

    Example:
        >>> input_variants = ['chr6-161006172-T-G']
        >>> from_genome = "hg19"
        >>> dest_genome = "hg38"
        >>> lifted_variants = lift_over_variants(input_variants, from_genome, dest_genome)
        >>> print(lifted_variants)
        ['chr6-160585140-T-G']
    """
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json",
        "User-Agent": user_agent,
    }

    # Parse input strings into DnaChange objects
    dna_changes = [DnaChange(*variant.split("-")) for variant in variants]

    # Prepare the payload
    variants = [
        {
            "chr": dna_change.chr,
            "pos": dna_change.pos,
            "ref": dna_change.ref,
            "alt": dna_change.alt,
        }
        for dna_change in dna_changes
    ]

    print(f"this is my payload: {json.dumps(variants)}")

    if (username is not None) and (api_key is not None):
        auth = (username, api_key)
        if use_netrc:
            _save_credentials_to_netrc(endpoint_url, username, api_key)
    elif use_netrc:
        username, _, api_key = _read_netrc_credentials(endpoint_url)
        if (username is not None) and (api_key is not None):
            auth = (username, api_key)
        else:
            auth = None
    else:
        auth = None

    params = {
        "from": from_genome,
        "dest": dest_genome,
    }

    batch_size = 500
    # Make the POST request

    result = []
    for i in tqdm(range(0, len(variants), batch_size)):
        # Prepare data for API request
        chunk = variants[i : i + batch_size]
        response = requests.post(
            endpoint_url, params=params, json=chunk, headers=headers, auth=auth
        )

        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            # Assuming the response contains a list of lifted variants
            lifted_variants = response.json()
            result.extend(
                (
                    f"{element.get('chr', '')}-{element.get('pos', '')}-{element.get('ref', '')}-{element.get('alt', '')}"
                    if element.get("chr")
                    else ""
                )
                for element in lifted_variants
            )
        else:
            # Print an error message if the request was not successful
            print(f"Error: {response.status_code}, {response.text}")
            return []
    return result


def parse_hgvs(
    hgvs: List[str],
    batch_size: int = 500,
    username: str = None,
    api_key: str = None,
    use_netrc: bool = True,
    endpoint_url: str = "https://api.genebe.net/cloud/api-public/v1/hgvs",
) -> List[str]:
    """
    Parses a list of genetic variants encoded in HGVS.

    Args:
        hgvs (List[str]): A list of genetic variants in hgvs format.
            Supports .n:, .c:, .g: and .m:. Look at examples below.
        batch_size (int, optional): The size of each batch for processing variants.
            Defaults to 500. Must be smaller or equal 1000.
        username (str, optional): The username for authentication.
            Defaults to None.
        api_key (str, optional): The API key for authentication.
            Defaults to None.
        use_netrc (bool, optional): Whether to use credentials from the user's
            .netrc file for authentication. Defaults to True.
        endpoint_url (str, optional): The API endpoint for parsing hgvs.
            Defaults to 'https://api.genebe.net/cloud/api-public/v1/hgvs'.

    Returns:
        List[str]: A list of dictionaries containing annotation information
        for each variant. Check the current documentation
        on https://genebe.net/about/api

    Example:
        >>> hgvs_variants = ["'NM_000277.2:c.1A>G"]
        >>> variants = parse_hgvs(hgvs_variants)
        >>> print(variants)
        ['12-102917129-AT-AC', '12-102852850-GA-G']

    Note:
        - The number of the elements in returned list is always equal to the number of queries.
        If some position cannot be parsed, an empty string is returned.
    """

    # logging in
    auth = None
    result = []

    for i in tqdm(range(0, len(hgvs), batch_size)):
        # Prepare data for API request
        chunk = hgvs[i : i + batch_size]

        logging.debug("Querying for " + json.dumps(chunk))
        # Make API request
        response = requests.post(
            endpoint_url,
            json=chunk,
            headers={
                "Accept": "application/json",
                "Content-Type": "application/json",
                "User-Agent": user_agent,
            },
            auth=auth,
        )

        # Check if request was successful
        if response.status_code == 200:
            api_results_raw = response.json()
            logging.debug("Api reqult raw" + json.dumps(api_results_raw))
            api_results = [
                (
                    f"{element.get('chr', '')}-{element.get('pos', '')}-{element.get('ref', '')}-{element.get('alt', '')}"
                    if "chr" in element and "pos" in element
                    else ""
                )
                for element in api_results_raw
            ]

            logging.debug("Backend result is " + json.dumps(api_results))

            # Append API results to annotated_data
            result.extend(api_results)
        elif response.status_code == 429:
            _handle_too_many_requests()
        else:
            logging.error(
                f"Got response with code {response.status_code} with body "
                + response.text
            )
            raise WrongServerResponseError(
                f"Got response with code {response.status_code} with body "
                + response.text
            )
    return result


def annotate_dataframe_variants(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    """
    Annotates genetic variants in a DataFrame using the 'annotate_variants_list_to_dataframe' function.

    Args:
        df (pd.DataFrame): Input DataFrame with genetic variant information. Must contain columns
        ["chr", "pos", "ref", "alt"]
        **kwargs: Additional keyword arguments to pass to 'annotate_variants_list_to_dataframe'.

    Returns:
        pd.DataFrame: A new DataFrame with additional annotation information.

    Example:
        >>> df = pd.DataFrame({'chr': ['6', '22'], 'pos': [160585140, 28695868], 'ref': ['T', 'AG'], 'alt': ['G', 'A']})
        >>> annotated_df = annotate_dataframe_variants(df, genome='hg38',use_ensembl=False,use_refseq=True, genome='hg38', flatten_consequences=True)
        >>> print(annotated_df)
    """
    # Ensure required columns are present in the DataFrame
    required_columns = ["chr", "pos", "ref", "alt"]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns in DataFrame: {missing_columns}")

    # Create a list of variant strings by concatenating columns
    variant_strings = [
        f"{row['chr']}-{row['pos']}-{row['ref']}-{row['alt']}"
        for _, row in df.iterrows()
    ]

    # Annotate variants using 'annotate_variants_list_to_dataframe'
    annotation_df = annotate_variants_list_to_dataframe(
        variants=variant_strings, **kwargs
    )
    annotation_df = annotation_df.drop(columns=["chr", "pos", "ref", "alt"])

    # Join the original DataFrame with the annotation results
    result_df = pd.concat([df.reset_index(drop=True), annotation_df], axis=1)

    return result_df


def whoami(
    username: str = None,
    api_key: str = None,
    use_netrc: bool = True,
    endpoint_url: str = "https://api.genebe.net/cloud/api-public/v1/whoami",
) -> List[str]:
    """
    Information about me: usage statistics, limits etc.

    Args:
        username (str, optional): The username for authentication.
            Defaults to None.
        api_key (str, optional): The API key for authentication.
            Defaults to None.
        use_netrc (bool, optional): Whether to use credentials from the user's
            .netrc file for authentication. Defaults to True.
        endpoint_url (str, optional): The API endpoint to use.
            Defaults to 'https://api.genebe.net/cloud/api-public/v1/whoami'.

    Returns:
        An object. Check the current documentation
        on https://genebe.net/about/api

    """

    # logging in
    if (username is not None) and (api_key is not None):
        auth = (username, api_key)
        if use_netrc:
            _save_credentials_to_netrc(endpoint_url, username, api_key)
    elif use_netrc:
        username, _, api_key = _read_netrc_credentials(endpoint_url)
        if (username is not None) and (api_key is not None):
            auth = (username, api_key)
        else:
            auth = None
    else:
        auth = None

    # Make API request
    response = requests.get(
        endpoint_url,
        headers={
            "Accept": "application/json",
            "Content-Type": "application/json",
            "User-Agent": user_agent,
        },
        auth=auth,
    )

    # Check if request was successful
    if response.status_code == 200:
        data = response.json()
        return data
    else:
        logging.warning("Error:", response.status_code, response.text)
        return {}
