"""VCF variant site data storage class.
"""

from typing import Dict, List, Union


class VariantSite:
    """Extracts and stores data about a VCF variant site.

    Attributes:
        chrom:
            String chromosome identifier from first column of the VCF row.
        pos:
            Integer chromosomal position from second column of the VCF row.
        ids:
            String list of site identifiers from third column of the VCF row.
            This will be '.' if no identifiers were assigned to the site.
        ref:
            String reference allele from fourth column of the VCF row.
        alts:
            List of string alternate alleles from fifth column of the VCF row.
        genotype_calls:
            Dictionary of values from the genotype fields of the VCF rows.
            The keys are based on the ninth (FORMAT) column.
        depth:
            Integer total depth at site combined across samples of the VCF row.
            This is retrieved from the seventh (INFO) column if available, or
            summed based on genotype fields if not.
        sample_names:
            String sample name list from the genotype fields of the VCF rows.
    """

    def __init__(self, vcf_line: List[str], sample_names: List[str]):
        self.chrom = vcf_line[0]
        self.pos = int(vcf_line[1])
        self.ids = vcf_line[2].split(';')
        self.ref = vcf_line[3]
        self.alts = vcf_line[4].split(',')
        genotype_format = vcf_line[8].split(':')
        self.genotype_calls: Dict[str, Dict] = {
            sample: self._parse_genotype_field(call, genotype_format)
            for sample, call in zip(sample_names, vcf_line[9:])
        }
        self.depth = self._retrieve_depth(vcf_line[7])

    def _retrieve_depth(self, info: str):
        depth = 0
        parsed_info = self._parse_info(info)
        if 'DP' in parsed_info:
            depth = parsed_info['DP']
        else:
            for call in self.genotype_calls.values():
                if 'DP' in call and call['DP'] != '.':
                    depth += int(call['DP'])
        return depth

    def _parse_genotype_field(self, call: str, genotype_format: List[str]):
        """Parse a genotype field based on processed FORMAT string.

        Args:
            call:
                String containing colon-separated values from a VCF genotype
                field for a single site and sample.
            genotype_format:
                List of types of data (GT, GQ, DP...) included in VCF genotype
                fields for a site. The list should match the contents and order
                of the (colon-separated) contents of the call argument.

        Returns:
            Dictionary of data subfields from a genotype field.
            The colon-separated elements of the FORMAT field serve as keys,
            while the colon-separated elements of the genotype field are values.
        """
        parsed_call: Dict[str, Union[str, int]] = {
            field: value
            for field, value in zip(genotype_format, call.split(':'))
        }
        return parsed_call

    def _parse_info(self, info: str):
        """Parse the additional information (INFO) field of a VCF site.

        Args:
            String containing INFO field (semicolon separated keys with optional
            values in a key=data format).

        Returns:
            Dictionary of information from an INFO field.
        """
        parsed_info = {}
        for subfield in info.split(';'):
            key_value = subfield.split('=')
            if len(key_value) == 1:
                parsed_info[key_value[0]] = True
            else:
                parsed_info[key_value[0]] = key_value[1]
        if 'DP' in parsed_info:
            parsed_info['DP'] = int(parsed_info['DP'])
        return parsed_info

    @property
    def sample_names(self):
        """Retrieve sample names.

        These are IDs from the source VCF header.
        """
        return list(self.genotype_calls.keys())

    def has_only_snvs(self):
        """Check if the site contains exclusively single nucleotide variants.

        The site is assumed to be SNV-only when both the reference allele and
        every alternate allele are one character long (indels would be longer).
        """
        return len(self.ref) != 1 or any(len(alt) != 1 for alt in self.alts)
