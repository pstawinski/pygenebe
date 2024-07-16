## GeneBe Utils
PyGeneBe: A Python client seamlessly integrating with the GeneBe platform, offering efficient annotation of genetic variants through its API, while supporting pandas, VCF file formats, and HGVS parsing

Using this client, you can easily annotate your DNA variants with the GeneBe API. Annotations include:
* Gene, transcript, and effect
* ClinVar phenotype
* GnomAD frequency
* ACMG score
* ... if you need more, please let me know

### Usage

For more information about the usage, go to the https://pygenebe.readthedocs.io/en/latest/ documentation.

#### Command line usage

Check current options using `--help` switch

```
genebe --help
genebe annotate --help
```


##### Annotating VCF using `annotate` command
GeneBe client allows you to annotate your VCF file with ease. Use the following command:

```
genebe annotate --input input.vcf.gz --output output.vcf.gz
```

Remember that your VCF file must be in a single allelic format! Utilize bcftools (https://samtools.github.io/bcftools/) to split the file. The output VCF will contain additional fields.

To use VCF annoation you have to have `cyvcf2` package installed. Take a look at the Installation section below.

If your VCF file is large (over 10.000 variants), you may encounter request limits. To avoid this, create a GeneBe account with an API Key and provide your login/key using the --username and --api-key arguments. You can always check your limits with the account command. Update your annotation command as follows:

```
genebe annotate --input input.vcf.gz --output output.vcf.gz --username your_username --api-key your_api_key
```

For more information call

```
genebe annotate --help
```

##### Using `account` command

The account command displays information about your request history statistics and limits. To check your limits without specifying a username and API key, run:

```
genebe account
```

Alternatively, if you have a GeneBe account with an API key, use the following command:

```
genebe account --username your_username --api-key your_api_key
```

Replace "your_username" and "your_api_key" with your GeneBe account credentials.

For more details and options, you can refer to the help documentation:

```
genebe account --help
```



#### Python usage

GeneBe makes annotating DNA variants in pandas dataframe easy.

```python
import genebe as gnb
import pandas as pd

input_variants = ['7-69599651-A-G']

# annotate variants with transcripts etc. output as a list; use .netrc for user login and api key
list = gnb.annotate(input_variants,
use_ensembl=True,
    use_refseq=False,
    genome="hg38",
    batch_size=500,
    use_netrc=True,
    output_format="list")


# output as a pandas dataframe, flat
df = gnb.annotate(input_variants,
    use_ensembl=True,
    use_refseq=False,
    genome="hg38",
    flatten_consequences=True,
    batch_size=500,
    use_netrc=True,
    output_format="dataframe")


# parse HGVS, SPDI or other
input_variants_parse = [
    "chrX:153803771:1:A",
    "22 28695868 AG A",
    "22-28695869--G",
    "22-28695869-G-",
    "NM_000277.2:c.1A>G",
    "NM_000277.2:c.2T>C",
    "AGT M259T",
    "rs1228544607"]
parsed_variants = gnb.parse_variants(input_variants_parse, genome="hg38")

# annotate existing dataframe, using it's chr, pos, ref, alt columns and adding new columns
df = pd.DataFrame({'chr': ['6', '22'], 'pos': [160585140, 28695868], 'ref': ['T', 'AG'], 'alt': ['G', 'A']})
annotated_df = gnb.annotate(df,
    genome='hg38',
    use_ensembl=False,
    use_refseq=True,
    flatten_consequences=True,
    output_format="dataframe")


# lift over variants from hg19 to hg38
input_variants = ['chr6-161006172-T-G']
from_genome = "hg19"
dest_genome = "hg38"
lifted_variants = gnb.lift_over_variants(input_variants, from_genome, dest_genome)


```

If you want to annotate thousands of variants, please log in to https://genebe.net, generate an API Key, and provide it using `username` and `api_key` or using the `.netrc` file.

Find out more usage examples in the `examples` directory.

### Installation
You can install GeneBe Utils using pip:

```
pip install genebe
```

If you wish to install faster `mmh3` implementation or use the option of annotating vcf files install using:

```
pip install genebe[cpp]
```

or install modules

```
pip install cyvcf2
pip install mmh3
```

in the environment.

This step will require build tools installed on your computer.

### Docker
There is a dockerized version of this package, available at https://hub.docker.com/r/genebe/pygenebe .

Usage example, reading from file `input.vcf` and writing output to `stdout`:
```
docker run -v input.vcf:/tmp/input.vcf --rm genebe/pygenebe:0.0.14 genebe annotate --input /tmp/input.vcf --output /dev/stdout
```

### Limits
If you wish to annotate thousands of variants, please log in to https://genebe.net, generate an API Key, and provide it using username and api_key.

The number of daily requests from a single IP is restricted to prevent abuse and excessive resource consumption on our server. Account holders with an API Key enjoy significantly higher limits (in the tens of thousands). If you require a higher daily request limit, please reach out to us via the https://genebe.net .

### Troubleshooting and issues
Experiencing issues? Follow these steps:

1. Check Existing Issues:

* If you encounter problems, explore existing issues on GitHub https://github.com/pstawinski/pygenebe for possible solutions.

2. Report New Issues:

* Unable to find a resolution? You can report the problem by creating a new issue with a clear description and details on https://github.com/pstawinski/pygenebe.

Your feedback is crucial for improving GeneBe client. Thank you for contributing to the community!

### Other

For more information about GeneBe, visit GeneBe website, https://genebe.net .




