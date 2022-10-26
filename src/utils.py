"""Utility functions for vcf2loc.
"""

import re
from typing import List, Union


def natural_sort_key(input_string: str) -> List[Union[int, str]]:
    """Key function for natural sorting of strings.

    Natural sorting treats multi-digit integers as a single value to be ordered
    numerically. This results in orderings like ['s1', 's2', 's10'] instead of
    the strictly "ASCIIbethical" ['s1', 's10', 's2'].

    Args:
        input_string:
            One of the alphanumeric strings to be sorted.

    Returns:
        List of (sub)strings and integers from the original string, to be used
        as sorting keys in the order they appear in the list.
    """
    return [
        int(substring) if substring.isdigit() else substring
        for substring in re.split(r'(\d+)', input_string)
    ]
