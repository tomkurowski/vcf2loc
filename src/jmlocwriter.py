"""Class used for writing JoinMap locus genotype (loc) files.
"""

import re
import sys
from typing import List

from src.jmmarker import JmMarker

SUPPORTED_POPULATION_TYPES = ['CP']


class JmLocWriter:
    """Opens a JoinMap locus genotype file and writes to it.

    Usable as a context manager (within with-statements).
    Markers can be written one-by-one by calling the write_marker method.
    The header metadata at the start of the file and the individual name list
    at the end of the file are handled automatically.

    Attributes:
        output_loc_filepath:
            String path (including file name) to the output JoinMap loc file.
        population_name:
            String name of the population. This can be up to 20 characters long
            and cannot contain spaces.
        population_type:
            String JoinMap population type ('CP').
        individual_names:
            List of strings containing names of individuals. These can be up
            to 20 characters long and cannot contain whitespace characters.
            May include parent names, which will be ignored.
        parent_a_name:
            String name of the first parent.
        parent_b_name:
            String name of the second parent.
        marker_count:
            Number of markers written so far.
        output_loc:
            Output JoinMap loc file handle.
    """

    def __init__(
        self,
        output_loc_filepath: str,
        population_name: str,
        population_type: str,
        parent_a_name: str,
        parent_b_name: str,
        individual_names: List[str]
    ):
        self.output_loc_filepath = output_loc_filepath
        self.population_name = population_name
        self.population_type = population_type
        self.individual_names = individual_names.copy()
        self.parent_a_name = parent_a_name
        self.parent_b_name = parent_b_name
        self.marker_count = 0
        self._validate_population_type(population_type)
        self.output_loc = open(self.output_loc_filepath, 'w', encoding='utf-8')
        self._write_loc_header()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        """Finish writing to JoinMap loc file."""
        self._write_individuals()
        self.output_loc.close()
        self._finalize_marker_count()

    def write_marker(self, marker: JmMarker):
        """Write a JoinMap marker row into the output file.

        Args:
            marker: JmMarker object.
        """
        if marker.segregation_type is None:
            print(marker.marker_name, file=self.output_loc)
        else:
            print(
                ' '.join([marker.marker_name, marker.segregation_type]),
                file=self.output_loc
            )
        ordered_gts = [
            marker.genotype_codes[name] for name in self.individual_names
            if name not in [self.parent_a_name, self.parent_b_name]
        ]
        print('   ' + ' '.join(ordered_gts), file=self.output_loc)
        self.marker_count += 1

    def _write_loc_header(self):
        """Write a JoinMap loc file header."""
        self._validate_population_name(self.population_name)
        child_count = len(self.individual_names)
        if self.parent_a_name in self.individual_names:
            child_count -= 1
        if self.parent_b_name in self.individual_names:
            child_count -= 1
        print(f'name = {self.population_name}', file=self.output_loc)
        print(f'popt = {self.population_type}', file=self.output_loc)
        print('nloc = MARKER_COUNT_PLACEHOLDER', file=self.output_loc)
        print(f'nind = {child_count}', file=self.output_loc)

    def _write_individuals(self):
        """Write a list of individual (sample) names to a JoinMap loc file."""
        print('individual names:', file=self.output_loc)
        for name in self.individual_names:
            if name not in [self.parent_a_name, self.parent_b_name]:
                self._validate_individual_name(name)
                print(name, file=self.output_loc)

    def _finalize_marker_count(self):
        """Replace the loc file header placeholder with the final marker count.
        """
        output = None
        with open(self.output_loc_filepath, 'r', encoding='utf-8') as outloc:
            output = outloc.read()
        output = output.replace(
            'MARKER_COUNT_PLACEHOLDER', str(self.marker_count)
        )
        with open(self.output_loc_filepath, 'w', encoding='utf-8') as outloc:
            outloc.write(output)

    def _validate_individual_name(self, name: str):
        if len(name) > 20:
            sys.exit(
                "Error: Individual names cannot be longer than 20 characters."
            )
        if re.search(r'\s', name):
            sys.exit(
                "Error: Individual names cannot contain whitespace characters."
            )

    def _validate_population_name(self, name: str):
        if len(name) > 20:
            sys.exit(
                "Error: Population name cannot be longer than 20 characters."
            )
        if re.search(r'\s', name):
            sys.exit(
                "Error: Population name cannot contain whitespace characters."
            )

    def _validate_population_type(self, population_type: str):
        if population_type not in SUPPORTED_POPULATION_TYPES:
            sys.exit(f"Error: Unsupported population type '{population_type}'.")
