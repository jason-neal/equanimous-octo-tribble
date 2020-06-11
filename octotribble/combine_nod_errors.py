import os
from glob import glob
from itertools import product

import numpy as np
from astropy.io import fits


def combine_errors(chips=range(1, 5), save_path=None):
    """Assumes the current directory contains the files for the nod cycle."""
    # Update header
    from datetime import datetime
    for chip in chips:
        # print("Chip", chip)
        ms_files = sorted(glob("C*_{0}.*ms.fits.*".format(chip)))
        norm_files = sorted(glob("C*_{0}.*ms.norm.fits.*".format(chip)))
        norm_sum_file = sorted(glob("C*_{0}.*ms.norm.sum.fits.*".format(chip)))

        assert len(ms_files) in (6, 8)
        assert len(norm_files) in (6, 8)
        assert len(norm_sum_file) == 1
        norm_sum_file = norm_sum_file[0]
        # Check ms_file and norm_file match the start of norm_sum_file
        assert any(os.path.basename(norm_sum_file)[0:29] in name for name in ms_files)
        assert any(os.path.basename(norm_sum_file)[0:29] in name for name in norm_files)

        # Combine Spectrum errors quadratically
        square_error = np.zeros(1024)
        for ii, (mf, nf) in enumerate(zip(ms_files, norm_files)):
            # Get relevant spectra from the fits files.
            msdata = fits.getdata(mf)
            normdata = fits.getdata(nf)
            assert msdata.shape == (3, 1, 1024)
            assert normdata.shape == (3, 1, 1024)

            ms_opt = msdata[0, 0, :].squeeze()
            norm_opt = normdata[0, 0, :].squeeze()
            ms_error = msdata[2, 0, :].squeeze()

            # Normalize the Error Spectrum
            norm_error = ms_error * norm_opt / ms_opt

            square_error += (norm_error * norm_error)

        # Squareroot and divide by number of spectra
        error_spectrum = np.sqrt(square_error) / (ii + 1)

        sumdata, header = fits.getdata(norm_sum_file, header=True)
        sum_spectrum = sumdata[0, 0, :].squeeze()

        assert error_spectrum.shape == sum_spectrum.shape

        header["history"] = "Errors combined on {}".format(datetime.now())
        header["comment"] = "Errors normalized, added in quadrature, then divide by number of nods"

        # Save sum_spectrum and error_spectrum together
        save_filename = "sum.witherr".join(norm_sum_file.split("sum"))
        print("Save filename = ", save_filename)

        if save_path is not None:
            # Change save location
            save_filename = os.path.join(save_path, os.path.basename(save_filename))

        save_errors(save_filename, spectrum=sum_spectrum, error=error_spectrum, header=header)


def save_errors(filename, spectrum, error, header):
    new_hdul = fits.HDUList()
    tbhdu = fits.BinTableHDU.from_columns(
        [fits.Column(name="Extracted_OPT", format="1D", array=spectrum),
         fits.Column(name="Error_OPT", format="1D", array=error)])
    tbhdu.header.extend(header)
    new_hdul.append(tbhdu)

    new_hdul.writeto(filename, output_verify="ignore", overwrite=True)
    return None


if __name__ == "__main__":
    orders = [45, 48]
    obsnames = ["CAL-OB1", "CAL-OB2", "CAL-OB3", "OB1", "OB2", "OB3"]
    chips = range(1, 5)
    combo = product(orders, obsnames)

    save_path = os.getcwd()

    for order, obsname in combo:
        print(order, obsname)

        path = f"/home/jneal/Phd/Collaborations/eta_tel_Hagelberg/data/good_reductions/order{order}/{obsname}/"
        os.chdir(path)

        combine_errors(save_path=save_path)

    os.chdir(save_path)
