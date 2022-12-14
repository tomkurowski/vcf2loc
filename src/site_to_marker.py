"""Functions for generating a JmMarker object based on a VariantSite object.
"""

from typing import Dict, List, Union

from src.variantsite import VariantSite
from src.jmmarker import JmMarker

_HOM_GENOTYPES = {'0/0', '1/1'}
_HET_GENOTYPES = {'0/1', '1/0'}

def is_potential_marker(
    site: VariantSite,
    population_type: str,
    parent_a: str,
    parent_b: str
):
    """Check if the site can serve as a marker based on parental variation.

    To serve as a useful marker, each of the parents needs to have one of a
    set of known genotype calls at the site.

    For the F2 and RIx population types, if both parents have a genotype call at
    the site, both calls have to be homozygous for different alleles. If only
    parent A has a genotype call, it has to be homozygous and at least one
    genotype call in the population has to be different from the parent A call.

    For the CP population type, at least one of the parents must have
    a heterozygous genotype, while the other can have either a homozygous or
    a heterozygous genotype.

    Args:
        site:
            VariantSite object generated based on a row of a VCF file.
        population_type:
            String JoinMap population type.
        parent_a:
            String name of first parent.
        parent_b:
            String name of second parent.
    """
    parent_a_gt = site.genotype_calls[parent_a]['GT']
    if population_type == 'CP':
        parent_b_gt = site.genotype_calls[parent_b]['GT']
        return (
            (parent_a_gt in _HET_GENOTYPES and parent_b_gt in _HET_GENOTYPES)
            or (parent_a_gt in _HOM_GENOTYPES and parent_b_gt in _HET_GENOTYPES)
            or (parent_a_gt in _HET_GENOTYPES and parent_b_gt in _HOM_GENOTYPES)
        )
    # If parent B was not provided, check if A is homozygous
    # and not identical to all other known samples.
    if parent_b is None:
        if parent_a_gt not in _HOM_GENOTYPES:
            return False
        gt_identical_to_parent_a = [
            site.genotype_calls[sample_name]['GT'] in { parent_a_gt, './.' }
            for sample_name in site.genotype_calls.keys()
        ]
        if all(gt_identical_to_parent_a):
            return False
        return True
    # Both parents were provided, check if they are homozygous and different
    parent_b_gt = site.genotype_calls[parent_b]['GT']
    return (
        parent_a_gt != parent_b_gt
        and parent_a_gt in _HOM_GENOTYPES
        and parent_b_gt in _HOM_GENOTYPES
    )


def site_to_marker(
    site: VariantSite,
    population_type: str,
    parent_a: str,
    parent_b: str,
    keep_invalid_calls: bool
) -> JmMarker:
    """Create and return a JmMarker object based on a VariantSite object.

    The first identifier from the site ID field is used as the marker name.
    If no identifiers are present (i.e. the ID field is '.'), the values of
    the CHROM and POS, separated by an underscore, are used for the marker name.

    Args:
        site:
            VariantSite object generated based on a row of a VCF file.
        population_type:
            String JoinMap population type.
        parent_a:
            String name of first parent.
        parent_b:
            String name of second parent.
        keep_invalid_calls:
            Boolean. If True, invalid / impossible genotype codes are recorded
            in the JmMarker object if the VCF calls in VariantSite support them.
            If False, the invalid codes are replaced by unknown codes.
    """
    genotypes = {
        sample_name: site.genotype_calls[sample_name]['GT']
        for sample_name in site.genotype_calls.keys()
    }
    if population_type == 'CP':
        segregation_type = _get_segregation_type(
            genotypes[parent_a], genotypes[parent_b]
        )
        genotype_codes = _convert_genotypes_cp(
            segregation_type,
            genotypes,
            parent_a,
            parent_b,
            site.sample_names,
            keep_invalid_calls
        )
        return JmMarker(
            marker_name=_generate_marker_name(site),
            genotype_codes=genotype_codes,
            population_type=population_type,
            segregation_type=segregation_type
        )
    # Generic approach for other population types.
    genotype_codes = _convert_genotypes(
        genotypes,
        parent_a,
        parent_b,
        site.sample_names
    )
    return JmMarker(
        marker_name=_generate_marker_name(site),
        genotype_codes=genotype_codes,
        population_type=population_type,
    )


def _generate_marker_name(site: VariantSite) -> str:
    if site.ids[0] != '.':
        return site.ids[0]
    return '_'.join([site.chrom, str(site.pos)])


def _get_segregation_type(genotype_parent_a: str, genotype_parent_b: str):
    if (
        genotype_parent_a in _HOM_GENOTYPES
        and genotype_parent_b in _HET_GENOTYPES
    ):
        segregation_type = '<nnxnp>'
    elif (
        genotype_parent_a in _HET_GENOTYPES
        and genotype_parent_b in _HOM_GENOTYPES
    ):
        segregation_type = '<lmxll>'
    elif (
        genotype_parent_a in _HET_GENOTYPES
        and genotype_parent_b in _HET_GENOTYPES
    ):
        segregation_type = '<hkxhk>'
    else:
        segregation_type = 'invalid'
    return segregation_type


def _convert_genotypes(
    genotypes: Dict[str, str],
    parent_a: str,
    parent_b: Union[str, None],
    sample_names: List[str]
):
    return {
        sample_name: _convert_genotype(
            genotypes[sample_name],
            genotypes[parent_a],
            None if parent_b is None else genotypes[parent_b]
        )
        for sample_name in sample_names
        if sample_name not in [parent_a, parent_b]
    }


def _convert_genotype(
    genotype: str,
    parent_a_genotype: str,
    parent_b_genotype: Union[str, None]
):
    genotype_code = '-'
    if genotype in _HOM_GENOTYPES:
        if genotype == parent_a_genotype:
            genotype_code = 'a'
        elif genotype == parent_b_genotype or parent_b_genotype is None:
            genotype_code = 'b'
    elif genotype in _HET_GENOTYPES:
        genotype_code = 'h'
    return genotype_code


def _convert_genotypes_cp(
    segregation_type: str,
    genotypes: Dict[str, str],
    parent_a: str,
    parent_b: str,
    sample_names: List[str],
    keep_invalid_calls: bool
):
    return {
        sample_name: _convert_genotype_cp(
            segregation_type,
            genotypes[sample_name],
            genotypes[parent_a],
            genotypes[parent_b],
            keep_invalid_calls
        )
        for sample_name in sample_names
        if sample_name not in [parent_a, parent_b]
    }


def _convert_genotype_cp(
    segregation_type: str,
    genotype: str,
    parent_a_genotype: str,
    parent_b_genotype: str,
    keep_invalid: bool
):
    genotype_code = '--'
    if genotype == './.':
        genotype_code = '--'
    if segregation_type == '<nnxnp>':
        genotype_code = _convert_nnxnp(
            genotype, parent_a_genotype, keep_invalid
        )
    elif segregation_type == '<lmxll>':
        genotype_code = _convert_lmxll(
            genotype, parent_b_genotype, keep_invalid
        )
    elif segregation_type == '<hkxhk>':
        genotype_code = _convert_hkxhk(genotype)
    return genotype_code


def _convert_nnxnp(genotype: str, parent_a_genotype: str, keep_invalid: bool):
    genotype_code = '--'
    if genotype in _HOM_GENOTYPES:
        if genotype == parent_a_genotype:
            genotype_code = 'nn'
        elif keep_invalid:
            # 'pp' is not a valid genotype for <nnxnp>
            genotype_code = 'pp'
        else:
            # Changing invalid genotype to unknown.
            genotype_code = '--'
    elif genotype in _HET_GENOTYPES:
        genotype_code = 'np'
    return genotype_code


def _convert_lmxll(genotype: str, parent_b_genotype: str, keep_invalid: bool):
    genotype_code = '--'
    if genotype in _HOM_GENOTYPES:
        if genotype == parent_b_genotype:
            genotype_code = 'll'
        elif keep_invalid:
            # 'mm' is not a valid genotype for <lmxll>
            genotype_code = 'mm'
        else:
            # Changing invalid genotype to unknown.
            genotype_code = '--'
    elif genotype in _HET_GENOTYPES:
        genotype_code = 'lm'
    return genotype_code


def _convert_hkxhk(genotype: str):
    genotype_code = '--'
    if genotype in _HOM_GENOTYPES:
        if genotype == '0/0':
            genotype_code = 'hh'
        elif genotype == '1/1':
            genotype_code = 'kk'
    elif genotype in _HET_GENOTYPES:
        genotype_code = 'hk'
    return genotype_code
