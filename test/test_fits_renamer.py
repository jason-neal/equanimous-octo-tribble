#test_fits_renamer.py

import pytest

from fits_renamer import main, char_replace


def test_char_replace():
    assert char_replace("test1:0.fits", ":","-") == "test1-0.fits"
    assert char_replace("test1:0.fits", ":","") == "test10.fits"


def test_main():
    # make file
    test_name = "test_make_file.tmp" 
    open(test_name, 'a').close()
    assert main(test_name,"_","-") == "test-make-file.tmp"
    # clean up
    os.remove("test-make-file.tmp")
    #subprocess.call(["rm", "test-make-file.tmp"])




if __name__ == "__main__":
	test_char_replace()
	test_main()