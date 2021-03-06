# test_fits_renamer.py
import os

import pytest

from octotribble.file_renamer import char_replace, main


def test_char_replace():
    assert char_replace("test1:0.fits", ":", "-") == "test1-0.fits"
    assert char_replace("test1:0.fits", ":", "") == "test10.fits"
    assert char_replace("test1_:_0.fits", "_", " ") == "test1 : 0.fits"
    assert char_replace("test1 0.fits", " ", "") == "test10.fits"


def test_main():
    test_name = "test_make_file.tmp"
    open(test_name, 'a').close()
    assert main(test_name, "_", "-") == "test-make-file.tmp"
    # clean up
    try:
        os.remove("test-make-file.tmp")
    except:
        os.remove("test_make_file.tmp")
    # subprocess.call(["rm", "test-make-file.tmp"])


if __name__ == "__main__":
    test_char_replace()
    test_main()
