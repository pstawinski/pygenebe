## GeneBe Utils
GeneBe: A Python client seamlessly integrating with the GeneBe platform, offering efficient annotation of genetic variants through its API, while supporting pandas, VCF file formats, and HGVS parsing

Using this client, you can easily annotate your DNA variants with the GeneBe API. Annotations include:
* Gene, transcript, and effect
* ClinVar phenotype
* GnomAD frequency
* ACMG score
* ... if you need more, please let me know

### Usage

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

input_variants = ['7-69599651-A-G']

# output as a list, with all transcripts
list = gnb.annotate_variants_list(input_variants,flatten_consequences = False)

# output as a pandas dataframe, flat
df = gnb.annotate_variants_list_to_dataframe(input_variants, flatten_consequences=True)

# parse HGVS
input_hgvs = ['NM_000277.2:c.1A>G']
parsed_variants = gnb.parse_hgvs(input_hgvs)

```

If you want to annotate thousands of variants, please log in to https://genebe.net, generate an API Key, and provide it using `username` and `api_key`.

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


### Other

For more information about GeneBe, visit GeneBe website, https://genebe.net .




