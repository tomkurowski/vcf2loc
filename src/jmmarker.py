"""JoinMap marker data storage class.
"""
from dataclasses import dataclass
from typing import Dict


@dataclass
class JmMarker:
    """Contains data for a JoinMap marker row.

    Attributes:
        marker_name:
            String name of the marker.
        genotype_codes:
            Dictionary of JoinMap genotype codes. The keys are individual names,
            and the values are genotype codes.
        segregation_type:
            String segregation type for CP population. These can be '<lmxll>',
            '<nnxnp>', or '<hkxhk>'. The '<abxcd>' and '<efxeg>' types are
            not used.
    """

    marker_name: str
    genotype_codes: Dict[str, str]
    segregation_type: str

    @property
    def known_count(self) -> int:
        """The number of known genotypes for this marker."""
        return len(self.genotype_codes) - self.unknown_count

    @property
    def unknown_count(self) -> int:
        """The number of unknown genotypes for this marker."""
        codes = list(self.genotype_codes.values())
        return codes.count('--')

    @property
    def unknown_fraction(self) -> float:
        """The fraction of unknown genotypes for this marker."""
        return self.unknown_count / len(self.genotype_codes)

    @property
    def homozygous_count(self) -> int:
        """The number of homozygous genotypes for this marker."""
        codes = list(self.genotype_codes.values())
        return (
            codes.count('hh')
            + codes.count('kk')
            + codes.count('ll')
            + codes.count('mm')
            + codes.count('nn')
            + codes.count('pp')
        )

    @property
    def homozygous_fraction(self) -> float:
        """The fraction of homozygous genotypes for this marker."""
        return self.homozygous_count / self.known_count

    @property
    def heterozygous_count(self) -> int:
        """The number of heterozygous genotypes for this marker."""
        codes = list(self.genotype_codes.values())
        return codes.count('hk') + codes.count('lm') + codes.count('np')

    @property
    def heterozygous_fraction(self) -> float:
        """The fraction of heterozygous genotypes for this marker."""
        return self.heterozygous_count / self.known_count
