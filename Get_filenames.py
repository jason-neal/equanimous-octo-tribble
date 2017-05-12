#!/usr/bin/env python
# -*- coding: utf8 -*-

"""Get Filenames of files that match regular expresion.

Possibly better to use the glob module.
"""
import os
import fnmatch
from typing import List


def get_filenames(path, regexp, regexp2=None):
    # type: (str, str, str) -> List[str]
    """Regexp must be a regular expression as a string.

    eg '*.ms.*', '*_2.*', '*.ms.norm.fits*'

    resexp2 is if want to match two expressions such as
    '*_1*' and '*.ms.fits*'
    """
    os.chdir(path)
    filelist = []
    for file in os.listdir('.'):
        if regexp2 is not None:  # Match two regular expresions
            if fnmatch.fnmatch(file, regexp) and fnmatch.fnmatch(file, regexp2):
                filelist.append(file)
        else:
            if fnmatch.fnmatch(file, regexp):
                filelist.append(file)
    filelist.sort()
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
