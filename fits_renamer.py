# Replace : with _ in fits files for windows

import argparse
import os
import subprocess

def _parser():
    """Take care of all the argparse stuff.
    :returns: the args
    """
    #parser = GooeyParser(description='Remove : from data files')
    parser = argparse.ArgumentParser(description='Rename fits files')
    parser.add_argument('filenames', type=list, help='List of Filenames')
    parser.add_argument("-x", "--cut", type=str, help="Character to Remove (cut)")
    parser.add_argument("-v", "--paste", type=str, help="Character to replace with (paste)")
    args = parser.parse_args()
    return args

def char_replace(fname, remove_char, replace_char=""):
    name_parts = fname.split(remove_char)
    new_name = replace_char.join(name_parts)
    return new_name


def main(filenames, cut=":", paste=""):
    for fname in filenames:
        new_name = char_replace(fname, cut, paste)
        #renaming file in shell
        os.rename(fname, new_name)
        #subprocess.call(["mv ",fname, new_name], )
        

if __name__ == '__main__':
    args = vars(_parser())
    #fname = args.pop('fname')
    opts = {k: args[k] for k in args}

    #main(fname, **opts)
    main(**opts)