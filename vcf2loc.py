#! /usr/bin/env python3
"""vcf2loc.py script for converting VCF files to JoinMap loc files.
"""
import argparse

from src import GbsVcfReader
from src import JmLocWriter
from src import site_to_marker, is_potential_marker
from src import natural_sort_key


parser = argparse.ArgumentParser(
    description="""Generate JoinMap loc file
based on vcf output from a GBS pipeline."""
)
parser.add_argument('input_vcf', type=str, help="input GBS vcf file")
parser.add_argument(
    '--output-loc', '-o',
    required=True,
    type=str,
    help="output JoinMap loc file"
)
parser.add_argument(
    '--name',
    required=False,
    default="GBS",
    type=str,
    help="output file population name; cannot contain spaces"
)
parser.add_argument(
    '--population-type', '-t',
    required=True,
    default=None,
    type=str,
    help="population type (supported: CP, F2, RIx)"
)
parser.add_argument(
    '--parent-a', '-a',
    required=True,
    default=None,
    type=str,
    help="first parent sample name"
)
parser.add_argument(
    '--parent-b', '-b',
    default=None,
    type=str,
    help="second parent sample name (optional except for CP population type)"
)
parser.add_argument(
    '--snvs',
    required=False,
    action='store_true',
    help="use only SNV sites"
)
parser.add_argument(
    '--min-dp',
    required=False,
    type=int,
    help="minimum DP (total depth, excluding parents) threshold for site"
)
parser.add_argument(
    '--u-threshold',
    required=False,
    default=1.0,
    type=float,
    help="maximum fraction of unknown calls per site"
)
parser.add_argument(
    '--het-threshold',
    required=False,
    default=1.0,
    type=float,
    help="maximum fraction (among known calls) of heterozygous calls per site"
)
parser.add_argument(
    '--hom-threshold',
    required=False,
    default=1.0,
    type=float,
    help="maximum fraction (among known calls) of homozygous calls per site"
)
parser.add_argument(
    '--keep-invalid-calls',
    required=False,
    action='store_true',
    default=False,
    help="keep calls which could not result from the parental genotypes"
)
parser.add_argument(
    '--natural-sort',
    required=False,
    action='store_true',
    default=False,
    help="apply a natural ordering to the individual names"
)

args = parser.parse_args()

if args.population_type == 'CP' and args.parent_b is None:
    parser.error("--parent-b / -b is required when using CP population type")

with GbsVcfReader(args.input_vcf) as invcf, \
     JmLocWriter(
        args.output_loc,
        args.name,
        args.population_type,
        args.parent_a,
        args.parent_b,
        sorted(invcf.sample_names, key=natural_sort_key)
        if args.natural_sort else invcf.sample_names
    ) as outloc:
    for site in invcf:
        if len(site.alts) > 1:
            # Skipping multiallelic sites.
            continue
        if args.snvs and site.has_only_snvs():
            # Keeping only SNV sites.
            continue
        if args.min_dp is not None:
            parent_b_dp = (
                0 if args.parent_b is None
                else site.genotype_calls[args.parent_a]['DP']
            )
            child_dp = (
                site.depth
                - site.genotype_calls[args.parent_a]['DP']
                - parent_b_dp
            )
            if child_dp < args.min_dp:
                # Skipping sites with low depth.
                continue
        if not is_potential_marker(
            site, args.population_type, args.parent_a, args.parent_b
        ):
            # Skipping sites which cannot serve as markers.
            continue
        marker = site_to_marker(
            site,
            args.population_type,
            args.parent_a,
            args.parent_b,
            args.keep_invalid_calls
        )
        if marker.unknown_fraction > args.u_threshold:
            # Skipping sites with too many unknown genotypes.
            continue
        if marker.homozygous_fraction > args.hom_threshold:
            # Skipping sites with too many homozygous genotypes.
            continue
        if marker.heterozygous_fraction > args.het_threshold:
            # Skipping sites with too many heterozygous genotypes.
            continue
        outloc.write_marker(marker)
