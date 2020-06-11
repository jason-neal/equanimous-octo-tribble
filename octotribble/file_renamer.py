#!/usr/bin/python
"""Motivation: Replace : with _ in fits files for use in windows.

# Author: Jason J. Neal <jason.neal@astro.up.pt>
"""
from __future__ import division, print_function

import argparse
import os


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    # parser = GooeyParser(description='Remove : from data files')
    parser = argparse.ArgumentParser(description='Rename fits files')
    parser.add_argument('filenames', nargs='+', type=str, help='List of Filenames')
    parser.add_argument("-x", "--cut", type=str, default=":", help="Character to Remove (cut)")
    parser.add_argument("-v", "--paste", type=str, default="", help="Character to replace with (paste)")
    return parser.parse_args()


def char_replace(fname, remove_char, replace_char=""):
    """Replace all of one character from a given string.

       inputs:
       fname: String of filename to replace characters from
       remove_char: Character to remove
    Keyword Arguments:

       - replace_char -- Character to replace with. Default = "" (blank)

    Examples:
    >>> char_replace("test_:.txt", ":", "1")
    'test_1.txt'

    >>> char_replace("test_adding_nothing.txt", "_")
    'testaddingnothing.txt'

    >>> char_replace("test_adding_space.txt", "_", " ")
    'testaddingnothing.txt'

    >>> char_replace("test removing space.txt", " ", "_" )
    'test_removing_space.txt'

    """
    name_parts = fname.split(remove_char)
    return replace_char.join(name_parts)


def main(filenames, cut=":", paste=""):
    """Replace characters in filenames.

    This script can be run from command line as
    python fits_renamer.py "test_filename.txt" -x "_" -v "-"

    # # TO DO: Could maybe exend this to take pairs to characters to replace.
    # # e.g. -x -, :, + -v _, _, -
    # # -x and -v would need to be same length

    """
    if isinstance(filenames, list):
        for fname in filenames:
            new_name = char_replace(fname, cut, paste)
            os.rename(fname, new_name)

    elif isinstance(filenames, str):
        new_name = char_replace(filenames, cut, paste)
        os.rename(filenames, new_name)
    return new_name


if __name__ == '__main__':
    args = vars(_parser())
    # fname = args.pop('fname')
    opts = {k: args[k] for k in args}

    # main(fname, **opts)
    main(**opts)
