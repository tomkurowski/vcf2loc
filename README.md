# vcf2loc

Convert VCF files generated by GBS pipelines like [GB-eaSy](https://github.com/dpwickland/GB-eaSy) or [TASSEL 5 GBSv2](https://bitbucket.org/tasseladmin/tassel-5-source/wiki/Tassel5GBSv2Pipeline) to the locus genotype file (loc-file) format used by [JoinMap](https://www.kyazma.nl/index.php/JoinMap/), which can be used to create linkage maps for experimental populations of diploid species.

Currently supports the CP, F2, and RIx population types. Multiallelic sites are ignored and phasing is not supported.

## Requirements

The script requires Python 3.7 or higher. No installation or dependencies outside the Python standard library are required.

## Usage

You can run the vcf2loc.py script file directly:

```sh
./vcf2loc.py -t CP -a PARENT_A -b PARENT_B -o gbs_output.loc gbs_input.vcf
```

Or (particularly on Windows):

```sh
python vcf2loc.py -t CP -a PARENT_A -b PARENT_B -o gbs_output.loc gbs_input.vcf
```

As in the above example, you need to provide the script with the type of population, the sample names of the two parental samples (these should match their names in the VCF file header), an output file name, and an (uncompressed) VCF input file.

You can view a help message, which details these and other (optional) arguments, through ```./vcf2loc.py --help``` or ```python vcf2loc.py --help```.

```
usage: vcf2loc.py [-h] --output-loc OUTPUT_LOC [--name NAME] --population-type POPULATION_TYPE --parent-a PARENT_A [--parent-b PARENT_B] [--snvs] [--min-dp MIN_DP] [--u-threshold U_THRESHOLD]
                  [--het-threshold HET_THRESHOLD] [--hom-threshold HOM_THRESHOLD] [--keep-invalid-calls] [--natural-sort]
                  input_vcf

Generate JoinMap loc file based on vcf output from a GBS pipeline.

positional arguments:
  input_vcf             input GBS vcf file

options:
  -h, --help            show this help message and exit
  --output-loc OUTPUT_LOC, -o OUTPUT_LOC
                        output JoinMap loc file
  --name NAME           output file population name; cannot contain spaces
  --population-type POPULATION_TYPE, -t POPULATION_TYPE
                        population type (supported: CP, F2, RIx)
  --parent-a PARENT_A, -a PARENT_A
                        first parent sample name
  --parent-b PARENT_B, -b PARENT_B
                        second parent sample name (optional except for CP population type)
  --snvs                use only SNV sites
  --min-dp MIN_DP       minimum DP (total depth, excluding parents) threshold for site
  --u-threshold U_THRESHOLD
                        maximum fraction of unknown calls per site
  --het-threshold HET_THRESHOLD
                        maximum fraction (among known calls) of heterozygous calls per site
  --hom-threshold HOM_THRESHOLD
                        maximum fraction (among known calls) of homozygous calls per site
  --keep-invalid-calls  keep calls which could not result from the parental genotypes
  --natural-sort        apply a natural ordering to the individual names
```

## Supported population type codes

Descriptions taken from the [JoinMap® 5 user manual by J.W. van Ooijen](https://www.kyazma.nl/docs/JM5Manual.pdf):

Type code|Description
---|---
```CP``` | a population resulting from a cross between two heterogeneously heterozygous and homozygous diploid parents, linkage phases originally possibly unknown
```F2``` | an F2 population: the result of selfing the F1 of a cross between two fully homozygous diploid parents
```RIx``` | a population of recombinant inbred lines in the x-th generation: the result of selfing an F2 with single seed descent; ```x``` must be specified: 2 <= ```x``` <= 99, ```RI2``` is equivalent to ```F2```

## Support & Contact information

Author & maintainer (Tomasz Kurowski): [t.j.kurowski@cranfield.ac.uk](mailto:t.j.kurowski@cranfield.ac.uk)
