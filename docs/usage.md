# Usage

## GeneBe Utils

GeneBe: A Python client seamlessly integrating with the GeneBe platform, offering efficient annotation of genetic variants through its API, while supporting pandas, VCF file formats, and HGVS parsing.

Using this client, you can easily annotate your DNA variants with the GeneBe API. Annotations include:

- Gene, transcript, and effect
- ClinVar phenotype
- GnomAD frequency
- ACMG score
- ... if you need more, please let us know

## Command Line Usage

Check current options using `--help` switch:

```bash
genebe --help
genebe annotate --help
```

### Annotating VCF using `annotate` command

GeneBe client allows you to annotate your VCF file with ease. Use the following command:

```bash
genebe annotate --input input.vcf.gz --output output.vcf.gz
```

Remember that your VCF file must be in a single allelic format! Utilize bcftools (https://samtools.github.io/bcftools/) to split the file. The output VCF will contain additional fields.

To use VCF annotation, you have to have the `cyvcf2` package installed. Take a look at the Installation section below.

If your VCF file is large (over 10,000 variants), you may encounter request limits. To avoid this, create a GeneBe account with an API Key and provide your login/key using the `--username` and `--api-key` arguments. You can always check your limits with the account command. Update your annotation command as follows:

```bash
genebe annotate --input input.vcf.gz --output output.vcf.gz --username your_username --api-key your_api_key
```


For more information, call:

```bash
genebe annotate --help
```

## Library usage

Example of annotating a single variant.
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

Find more usage examples on GitHub https://github.com/pstawinski/pygenebe/tree/main/examples
