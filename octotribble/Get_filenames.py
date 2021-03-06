#!/usr/bin/env python
# -*- coding: utf8 -*-

"""Get names of files that match regular expression.

Possibly better to use the glob module.
"""
import fnmatch
import os
from typing import List

# TODO: Try glob.glob


def get_filenames(path, regexp, regexp2=None, fullpath=False):
    # type: (str, str, str) -> List[str]
    """Regexp must be a regular expression as a string.

    eg '*.ms.*', '*_2.*', '*.ms.norm.fits*'

    regexp2 is if want to match two expressions such as
    '*_1*' and '*.ms.fits*'
    """
    current_path= os.getcwd()
    os.chdir(path)
    filelist = []
    for file in os.listdir('.'):
        if regexp2 is not None:  # Match two regular expressions
            if fnmatch.fnmatch(file, regexp) and fnmatch.fnmatch(file, regexp2):
                filelist.append(file)
        else:
            if fnmatch.fnmatch(file, regexp):
                filelist.append(file)
    filelist.sort()
    os.chdir(current_path)
    if fullpath:
        filelist = [os.path.join(path, f) for f in filelist]

    return filelist


def main():
    # type: () -> None
    """Some test examples."""
    path = "/home/jneal/data/BrownDwarfs-PedrosCode/HD30501-1/"

    list1 = get_filenames(path, "*.ms.*")
    for file in list1:
        pass  # print file

    list2 = get_filenames(path, "*.norm.*", "*_1.*")
    for file in list2:
        pass  # print file


if __name__ == '__main__':
    main()
