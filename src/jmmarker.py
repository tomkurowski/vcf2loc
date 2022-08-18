"""JoinMap marker data storage class.
"""
from dataclasses import dataclass
from typing import Dict, Optional


@dataclass
class JmMarker:
    """Contains data for a JoinMap marker row.

    Attributes:
        marker_name:
            String name of the marker.
        genotype_codes:
            Dictionary of JoinMap genotype codes. The keys are individual names,
            and the values are genotype codes.
        population_type:
            String JoinMap population type.
        segregation_type:
            String segregation type for CP population. These can be '<lmxll>',
            '<nnxnp>', or '<hkxhk>'. The '<abxcd>' and '<efxeg>' types are
            not used.
    """

    marker_name: str
    genotype_codes: Dict[str, str]
    population_type: str
    segregation_type: Optional[str] = None

    @property
    def known_count(self) -> int:
        """The number of known genotypes for this marker."""
        return len(self.genotype_codes) - self.unknown_count

    @property
    def unknown_count(self) -> int:
        """The number of unknown genotypes for this marker."""
        codes = list(self.genotype_codes.values())
        if self.population_type == 'CP':
            return codes.count('--')
        # Using '-' for all other population types ('u' or '.' not used).
        return codes.count('-')

    @property
    def unknown_fraction(self) -> float:
        """The fraction of unknown genotypes for this marker."""
        return self.unknown_count / len(self.genotype_codes)

    @property
    def homozygous_count(self) -> int:
        """The number of homozygous genotypes for this marker."""
        codes = list(self.genotype_codes.values())
        if self.population_type == 'CP':
            return (
                codes.count('hh')
                + codes.count('kk')
                + codes.count('ll')
                + codes.count('mm')
                + codes.count('nn')
                + codes.count('pp')
            )
        # 'a' and 'b' are the homozygous codes for non-CP diploid populations.
        return codes.count('a') + codes.count('b')

    @property
    def homozygous_fraction(self) -> float:
        """The fraction of homozygous genotypes for this marker."""
        return self.homozygous_count / self.known_count

    @property
    def heterozygous_count(self) -> int:
        """The number of heterozygous genotypes for this marker."""
        codes = list(self.genotype_codes.values())
        if self.population_type == 'CP':
            return codes.count('hk') + codes.count('lm') + codes.count('np')
        # 'h' is the heterozygous code for non-CP diploid populations.
        # 'c' and 'd' are not used as VCF files have no data on dominance.
        return codes.count('h')

    @property
    def heterozygous_fraction(self) -> float:
        """The fraction of heterozygous genotypes for this marker."""
        return self.heterozygous_count / self.known_count
