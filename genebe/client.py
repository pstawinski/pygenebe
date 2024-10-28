import pandas as pd
import requests
from typing import Dict, List, Union
import re
import json
import logging
import netrc
from tinynetrc import Netrc
from urllib.parse import urlparse
from .version import __version__
from pathlib import Path

try:
    from tqdm import tqdm
except ImportError:

    def tqdm(iterator, *args, **kwargs):
        return iterator


# Set up the User-Agent header
user_agent = f"GeneBe/{__version__} (requests)"
apikey_username_showed = False


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
        print(f"Warning: Invalid format in [{s}], not annotating.")
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
    machine_name = parsed_url.netloc.split(":")[0]
    return machine_name


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


def annotate(
    variants: Union[List[str], pd.DataFrame],
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
    output_format: str = "list",
    custom_annotations: List[str] = None,
) -> Union[List[Dict[str, object]], pd.DataFrame]:
    """
    Annotates genetic variants.

    Args:
        variants (Union[List[str], pd.DataFrame]): A list of genetic variants to be annotated
            or a pandas DataFrame containing the variants. If a list, the format should be
            chr-pos-ref-alt.
        genome (str): The genome version for annotation (e.g., 'hg38'). Defaults to 'hg38'.
        use_ensembl (bool): Whether to use Ensembl for annotation. Defaults to True.
        use_refseq (bool): Whether to use RefSeq for annotation. Defaults to True.
        flatten_consequences (bool): If set to False, return consequences as a list of objects.
            If set to True, only the most important consequence is returned in a flat form.
            Defaults to True.
        batch_size (int): The size of each batch for processing variants. Defaults to 500.
            Must be smaller or equal to 1000.
        username (str): The username for authentication. Defaults to None.
        api_key (str): The API key for authentication. Defaults to None.
        use_netrc (bool): Whether to use credentials from the user's .netrc file for authentication.
            Defaults to True.
        endpoint_url (str): The API endpoint for variant annotation.
            Defaults to 'https://api.genebe.net/cloud/api-public/v1/variants'.
        progress_bar (bool): Show progress bar. Defaults to True.
        omit_acmg (bool): Don't add ACMG scores in the output. Defaults to False.
        omit_csq (bool): Don't add consequences in the output. Defaults to False.
        omit_basic (bool): Don't add basic annotations (GnomAD frequencies etc) in the output. Defaults to False.
        omit_advanced (bool): Don't add advanced annotations (ClinVar etc) in the output. Defaults to False.
        omit_normalization (bool): Don't normalize variants. Use only if you are sure they are normalized already. Defaults to False.
        annotator (str): Which VEP implementation to use. Defaults to snpeff.
        custom_annotations (List[str]): A list of custom annotations to include in the output. Consult with the documentation for the list of available custom annotations.
        output_format (str): The desired format of the output, either 'list' for a list of
            dictionaries or 'dataframe' for a pandas DataFrame. Defaults to 'list'.
            If input is a dataframe, then ignored and returns a dataframe.

    Returns:
        Union[List[Dict[str, object]], pd.DataFrame]: If output_format is 'list', returns a list of
        dictionaries containing annotation information for each variant. If output_format is 'dataframe',
        returns a pandas DataFrame with the annotations.

    Raises:
        ValueError: If output_format is not one of 'list' or 'dataframe'.
        TypeError: If variants is neither a List[str] nor a pandas DataFrame.

    Example:
        >>> annotations = annotate(
        ...     ["7-69599651-A-G", "6-160585140-T-G"],
        ...     genome='hg38',
        ...     use_ensembl=True,
        ...     use_refseq=False,
        ...     batch_size=500,
        ...     username="user123@example.com",
        ...     api_key="apikey456",
        ...     use_netrc=False,
        ...     endpoint_url='https://api.genebe.net/cloud/api-public/v1/variants',
        ...     output_format="list"
        ... )
        >>> print(annotations)
        [{'chr': '7', 'pos': 69599651, 'ref': 'A', 'alt': 'G', 'annotation': '...'}, ...]

        >>> df = pd.DataFrame({"variant": ["7-69599651-A-G", "6-160585140-T-G"]})
        >>> annotations_df = annotate(
        ...     df,
        ...     genome='hg38',
        ...     use_ensembl=True,
        ...     use_refseq=False,
        ...     batch_size=500,
        ...     username="user123@example.com",
        ...     api_key="apikey456",
        ...     use_netrc=False,
        ...     endpoint_url='https://api.genebe.net/cloud/api-public/v1/variants',
        ...     output_format="dataframe"
        ... )
        >>> print(annotations_df)

    Note:
        - The number of elements in the returned list or rows in the DataFrame is always equal to the number of queries.
    """

    allowed_output_formats = ["list", "dataframe"]

    if output_format not in allowed_output_formats:
        raise ValueError(f"output_format must be one of {allowed_output_formats}")

    join_with_oryginal_df = False
    if isinstance(variants, pd.DataFrame):
        required_columns = ["chr", "pos", "ref", "alt"]
        missing_columns = [
            col for col in required_columns if col not in variants.columns
        ]
        if missing_columns:
            raise ValueError(
                f"Missing required columns in DataFrame: {missing_columns}"
            )

        # Store the original DataFrame for later use
        oryginal_df = variants
        join_with_oryginal_df = True

        # Create a list of variant strings by concatenating columns
        variant_strings = [
            f"{row['chr']}-{row['pos']}-{row['ref']}-{row['alt']}"
            for _, row in variants.iterrows()
        ]

        variants = variant_strings  # now we continue with the list
        output_format = "dataframe"

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
    global apikey_username_showed

    if api_key is None:
        if not apikey_username_showed:
            # Show it only once
            logging.warning(
                "You are not logged in to GeneBe. We recommend using username and api_key arguments. Find out more on https://genebe.net/api ."
            )
            apikey_username_showed = True
    else:
        if not apikey_username_showed:
            logging.info(f"I will try to log in as {username}")
            apikey_username_showed = True

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
        if custom_annotations:
            params["customAnnotations"] = ",".join(custom_annotations)

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

    if join_with_oryginal_df:
        annotation_df = pd.DataFrame(annotated_data)
        annotation_df = annotation_df.drop(columns=["chr", "pos", "ref", "alt"])

        # Join the original DataFrame with the annotation results
        result_df = pd.concat(
            [oryginal_df.reset_index(drop=True), annotation_df], axis=1
        )
        return result_df
    else:
        if output_format == "list":
            return annotated_data
        elif output_format == "dataframe":
            return pd.DataFrame(annotated_data)
        else:
            raise ValueError(f"output_format must be one of {allowed_output_formats}")


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
    DEPRECIATED. Use annotate function. Annotates a list of genetic variants.

    Args:
        variants (List[str]): A list of genetic variants to be annotated in the format chr-pos-ref-alt.
        genome (str, optional): The genome version for annotation (e.g., 'hg38'). Defaults to 'hg38'.
        use_ensembl (bool, optional): Whether to use Ensembl for annotation. Defaults to True.
        use_refseq (bool, optional): Whether to use RefSeq for annotation. Defaults to True.
        flatten_consequences (bool, optional): If set to False, return consequences as a list of objects.
            If set to True, only the most important consequence is returned in a flat form. Defaults to True.
        batch_size (int, optional): The size of each batch for processing variants. Defaults to 500.
            Must be smaller or equal to 1000.
        username (str, optional): The username for authentication. Defaults to None.
        api_key (str, optional): The API key for authentication. Defaults to None.
        use_netrc (bool, optional): Whether to use credentials from the user's .netrc file for authentication.
            Defaults to True.
        endpoint_url (str, optional): The API endpoint for variant annotation.
            Defaults to 'https://api.genebe.net/cloud/api-public/v1/variants'.
        progress_bar (bool, optional): Show progress bar. Defaults to True.
        omit_acmg (bool, optional): Don't add ACMG scores in the output. Defaults to False.
        omit_csq (bool, optional): Don't add consequences in the output. Defaults to False.
        omit_basic (bool, optional): Don't add basic annotations (GnomAD frequencies etc) in the output.
            Defaults to False.
        omit_advanced (bool, optional): Don't add advanced annotations (ClinVar etc) in the output.
            Defaults to False.
        omit_normalization (bool, optional): Don't normalize variants. Use only if you are sure they are normalized already.
            Defaults to False.
        annotator (str, optional): Which VEP implementation to use. Defaults to 'snpeff'.

    Returns:
        List[Dict[str, object]]: A list of dictionaries containing annotation information for each variant.

    Example:
        >>> variants = ["7-69599651-A-G", "6-160585140-T-G"]
        >>> annotations = annotate_variants_list(variants, use_ensembl=True,
        ...                                      use_refseq=False, genome='hg38',
        ...                                      batch_size=500, username="user123@example.com",
        ...                                      api_key="apikey456", use_netrc=False,
        ...                                      endpoint_url='https://api.genebe.net/cloud/api-public/v1/variants')
        >>> print(annotations)
        [{'chr': '7', 'pos': 69599651, ... }, {'chr': '6', 'pos': 160585140, ... }]

    Note:
        - The number of elements in the returned list is always equal to the number of queries.
    """

    # Call the annotate function with the variants and args
    annotated_variants = annotate(
        variants,
        genome=genome,
        use_ensembl=use_ensembl,
        use_refseq=use_refseq,
        flatten_consequences=flatten_consequences,
        batch_size=batch_size,
        username=username,
        api_key=api_key,
        use_netrc=use_netrc,
        endpoint_url=endpoint_url,
        progress_bar=progress_bar,
        omit_acmg=omit_acmg,
        omit_csq=omit_csq,
        omit_basic=omit_basic,
        omit_advanced=omit_advanced,
        omit_normalization=omit_normalization,
        annotator=annotator,
        output_format="list",
    )

    # Ensure the output is a pandas DataFrame
    if isinstance(annotated_variants, list):
        return annotated_variants
    else:
        raise ValueError("The output of annotate function is not a list.")


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
    """
    DEPRECIATED. Use annotate function.
    Annotates a list of variants and returns the annotations as a pandas DataFrame.

    Parameters:
        variants (List[str]): A list of variant identifiers.
        use_ensembl (bool, optional): Whether to use Ensembl for annotation. Defaults to True.
        use_refseq (bool, optional): Whether to use RefSeq for annotation. Defaults to True.
        genome (str, optional): The genome version for annotation (e.g., 'hg38'). Defaults to 'hg38'.
        batch_size (int, optional): The size of each batch for processing variants. Defaults to 100.
        flatten_consequences (bool, optional): If set to False, return consequences as a list of objects.
                                               If set to True, only the most important consequence is returned in a flat form. Defaults to True.
        username (str, optional): The username for authentication. Defaults to None.
        api_key (str, optional): The API key for authentication. Defaults to None.
        use_netrc (bool, optional): Whether to use credentials from the user's .netrc file for authentication. Defaults to True.
        endpoint_url (str, optional): The API endpoint for variant annotation. Defaults to 'https://api.genebe.net/cloud/api-public/v1/variants'.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the annotated variants.

    Raises:
        ValueError: If the output of the annotate function is not a pandas DataFrame.

    Example:
        >>> df = annotate_variants_list_to_dataframe(variants=["variant1", "variant2"])
        >>> print(df.head())
    """
    annotated_variants = annotate(
        variants,
        genome=genome,
        use_ensembl=use_ensembl,
        use_refseq=use_refseq,
        flatten_consequences=flatten_consequences,
        batch_size=batch_size,
        username=username,
        api_key=api_key,
        use_netrc=use_netrc,
        endpoint_url=endpoint_url,
        output_format="dataframe",
    )

    # Ensure the output is a pandas DataFrame
    if isinstance(annotated_variants, pd.DataFrame):
        return annotated_variants
    else:
        raise ValueError("The output of annotate function is not a pandas DataFrame.")


def lift_over_variants(
    variants: List[str],
    from_genome: str = "hg19",
    dest_genome: str = "hg38",
    username: str = None,
    api_key: str = None,
    use_netrc: bool = True,
    progress: bool = True,
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
        progress (bool, optional): Show progress bar. Defaults to True.
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

    vars = variants

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
    for i in tqdm(range(0, len(vars), batch_size)):
        # Prepare data for API request
        chunk = vars[i : i + batch_size]
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
                    if isinstance(element, dict) and element.get("chr")
                    else ""
                )
                for element in lifted_variants.get("variants")
            )
        else:
            # Print an error message if the request was not successful
            print(f"Error: {response.status_code}, {response.text}")
            return []
    return result


def parse_variants(
    ids: List[str],
    batch_size: int = 500,
    username: str = None,
    api_key: str = None,
    use_netrc: bool = True,
    progress: bool = True,
    genome: str = "hg38",
    endpoint_url: str = "https://api.genebe.net/cloud/api-public/v1/convert",
) -> List[str]:
    """
    Parses a list of genetic variants encoded in HGVS, SPDI, rs* or other format.

    Args:
        ids (List[str]): A list of genetic variants in hgvs format.
            Supports .n:, .c:, .g: and .m:. Look at examples below.
        batch_size (int, optional): The size of each batch for processing variants.
            Defaults to 500. Must be smaller or equal 1000.
        username (str, optional): The username for authentication.
            Defaults to None.
        api_key (str, optional): The API key for authentication.
            Defaults to None.
        use_netrc (bool, optional): Whether to use credentials from the user's
            .netrc file for authentication. Defaults to True.
        progress (bool, optional): Show progress bar. Defaults to True.
        endpoint_url (str, optional): The API endpoint for parsing hgvs.
            Defaults to 'https://api.genebe.net/cloud/api-public/v1/convert'.

    Returns:
        List[str]: A list of dictionaries containing annotation information
        for each variant. Check the current documentation
        on https://genebe.net/about/api

    Example:
        >>> input_variants = ["'NM_000277.2:c.1A>G", "22 28695868 AG A", "rs1228544607", "AGT Met259Thr"]
        >>> variants = parse_variants(input_variants)
        >>> print(variants)
        ['12-102917129-AT-AC', '12-102852850-GA-G', ']

    Note:
        - The number of the elements in returned list is always equal to the number of queries.
        If some position cannot be parsed, an empty string is returned.
    """

    # logging in
    auth = None
    result = []

    for i in tqdm(range(0, len(ids), batch_size)):
        # Prepare data for API request
        chunk = ids[i : i + batch_size]

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
            params={genome: genome},
            auth=auth,
        )

        # Check if request was successful
        if response.status_code == 200:
            api_results_raw = response.json()
            logging.debug("Api reqult raw" + json.dumps(api_results_raw))
            api_results = []

            for variants_dict in api_results_raw:
                if "variants" in variants_dict and variants_dict["variants"]:
                    variants = variants_dict["variants"]
                    variant = variants[0]
                    formatted_variant = f"{variant.get('chr', '')}-{variant.get('pos', '')}-{variant.get('ref', '')}-{variant.get('alt', '')}"
                    api_results.append(formatted_variant)
                else:
                    error_message = variants_dict.get(
                        "error", "No error message provided"
                    )
                    logging.warning(f"Warning: {error_message}")
                    api_results.append("")

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
    DEPRECIATED, use annotate

    Annotates genetic variants in a DataFrame using the 'annotate_variants_list_to_dataframe' function.

    Args:
        df (pd.DataFrame): Input DataFrame with genetic variant information. Must contain columns ["chr", "pos", "ref", "alt"].
        **kwargs: Additional keyword arguments to pass to 'annotate_variants_list_to_dataframe'.

    Returns:
        pd.DataFrame: A new DataFrame with additional annotation information.

    Example:
        >>> df = pd.DataFrame({'chr': ['6', '22'], 'pos': [160585140, 28695868], 'ref': ['T', 'AG'], 'alt': ['G', 'A']})
        >>> annotated_df = annotate_dataframe_variants(df, genome='hg38', use_ensembl=False, use_refseq=True, flatten_consequences=True)
        >>> print(annotated_df)

    Raises:
        ValueError: If the input DataFrame does not contain the required columns.
    """
    required_columns = {"chr", "pos", "ref", "alt"}
    if not required_columns.issubset(df.columns):
        raise ValueError(
            f"Input DataFrame must contain the following columns: {required_columns}"
        )

    # Prepare variants list in the required format
    variants = df.apply(
        lambda row: f"{row['chr']}-{row['pos']}-{row['ref']}-{row['alt']}", axis=1
    ).tolist()

    # Call the annotate_variants_list_to_dataframe function with the prepared variants and additional kwargs
    annotated_df = annotate_variants_list_to_dataframe(variants, **kwargs)

    return annotated_df


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
