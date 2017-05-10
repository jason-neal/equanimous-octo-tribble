"""Identify bad pixel in nod observations and interpolate over them."""
import logging
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt


def sigma_detect(nods, plot=True):
    """Detect the local pixels that are outside 4sigma from all nods.

    Local means 2 pixels either side.
    """
    if isinstance(nods, list):
        raise TypeError("Input an nod*pixel array please.")
    sig_clip = 4  # Sigma clipping Value.
    if nods.shape[0] > 8:
        raise ValueError("Too many nods (>8), check dimensions of input. ([nod, pixel])")
    # aviod mutation
    old_nods = np.empty_like(nods)
    old_nods[:] = np.nan   # set to nans to start iteration
    new_nods = np.empty_like(nods)
    new_nods[:] = nods

    bad_pixel_count = 0
    bad_pixel_record = []

    iteration = 0
    # Iterate untill no more bad pixels are replaced by nans during an iteration. or less than 5.
    while ((iteration < 5) & np.any(np.isnan(old_nods) != np.isnan(new_nods))):
        old_nods[:] = new_nods
        for pixel in range(new_nods.shape[1]):
            if (pixel < 2):
                near_pixels = new_nods[:, :5]      # First 5 pixels to do the 2 end pixels..
                grid_index = pixel
            elif pixel > (new_nods.shape[1] - 3):
                near_pixels = new_nods[:, -5:]     # Last 5 pixels to do the last 2 end pixels..
                grid_index = pixel - new_nods.shape[1]
            else:
                near_pixels = new_nods[:, slice(pixel - 2, pixel + 3)]
                grid_index = 2

            # ravel pixels near this picel output
            ravel_pixels = near_pixels.ravel()
            median = np.nanmedian(ravel_pixels)   # ignore nan values
            std = np.nanstd(ravel_pixels)         # ignore nan values

            # The values of this pixel for all nods, taken from near_pixels. Should be same as new_nods[:, pixel]
            this_pixel = near_pixels[:, grid_index]
            assert np.all((this_pixel == new_nods[:, pixel]) | np.isnan(this_pixel))

            # Find values outside the sig_clip level.
            sig = (this_pixel > (median + sig_clip * std)) | (this_pixel < (median - sig_clip * std))

            if np.any(sig):
                bad_nod = sig.nonzero()[0]
                if len(bad_nod) > 1:
                    logging.Warning("More then one nod has a bad pixel value here, in pixel #{}".format(pixel))
                for val in bad_nod:
                    bad_pixel_count += 1
                    bad_pixel_record += [(val, pixel, this_pixel[val])]
        for record in bad_pixel_record:
            new_nods[record[0], record[1]] = np.nan
        iteration += 1

    print("# Pixels outside {0}sigma = {1}".format(sig_clip, bad_pixel_count))
    # print("bad_pixel_record", bad_pixel_record)

    if plot:
        for i, nod in enumerate(nods):
            plt.plot(nod, label="nod {}".format(i))
        bad_pixel_x = [x[1] for x in bad_pixel_record]
        bad_pixel_y = [x[2] for x in bad_pixel_record]

        plt.plot(bad_pixel_x, bad_pixel_y, "x", label=">4 sigma")
        plt.legend()
        plt.title("Identifying bad pixels")
        plt.xlabel("pixel")
        plt.ylabel("norm flux")
        plt.show()

    return [pixel[0:2] for pixel in bad_pixel_record]


def interp_badpixels(nods, bad_pixels):
    """Linearly interpolate over nearby pixels.

    If it is at the end then just replace with the next pixel value.

    Parameters:
    nods: array
        A nod*pixel array of the pixel flux values.
    bad_pixels: list of lists of ints
        Index position [nod, pixel] of bad pixels to replace.

    Returns
    -------
    nods: array
        Array of nods with bad pixels interpolated over.
    """
    if isinstance(nods, list):
        raise TypeError("Input an nod*pixel array please.")

    if not isinstance(bad_pixels, list):
        raise TypeError("Bad_pixels must be a list.")

    if nods.shape[1] < 3:
        raise ValueError("The number of pixels to interpolate ({}) is very small."
                         "Are you sure you are doing the correct thing?".format(nods.shape[1]))

    # Warn about consecutive bad_pixels
    warn_consec_badpixels(bad_pixels, stop=False)
    output = np.empty_like(nods)
    output[:] = nods            # aviods mutation

    for pixel in bad_pixels:
        # Count any consecutive bad pixels
        bp_right = right_consec_search(pixel, bad_pixels)
        bp_left = left_consec_search(pixel, bad_pixels)

        if (bp_left + bp_right) > 5:
            logging.warning("Interpolating over more than 5 bad pixels in a row.")

        if (bp_left + bp_right + 1) == nods.shape[1]:
            replacement = nods[pixel[0], pixel[1]]           # Replace with self as all others are bad.
        elif pixel[1] == 0:                                  # first vaue
            replacement = nods[pixel[0], pixel[1] + 1 + bp_right]

        elif pixel[1] == (nods.shape[1] - 1):                # Last vaue
            replacement = nods[pixel[0], pixel[1] - 1 - bp_left]

        elif (nods.shape[1] - pixel[1] - 1) == bp_right:   # All bad pixels until the end
                    replacement = nods[pixel[0], pixel[1] - 1 - bp_left]

        elif pixel[1] == bp_left:                             # All bad pixels to the left
                    replacement = nods[pixel[0], pixel[1] + 1 + bp_right]
        else:
            x = range(0 - bp_left, 3 + bp_right)
            xp = [0 - bp_left, 2 + bp_right]
            fp = [nods[pixel[0], pixel[1] - 1 - bp_left], nods[pixel[0], pixel[1] + 1 + bp_right]]
            y = np.interp(x, xp, fp)
            replacement = y[1 + bp_left]

        output[pixel[0], pixel[1]] = replacement

    return output


def left_consec_search(pixel, bad_pixels):
    """Count number of consecutive bad pixels to the left of this pixel."""
    prev_pixel = [pixel[0], pixel[1] - 1]
    if (prev_pixel in bad_pixels) or (tuple(prev_pixel) in bad_pixels):
        # print("prev_pixel in recursion", prev_pixel)
        return left_consec_search(prev_pixel, bad_pixels) + 1
    else:
        return 0


def right_consec_search(pixel, bad_pixels):
    """Count number of consecutive bad pixels to the right of this."""
    next_pixel = [pixel[0], pixel[1] + 1]
    if (next_pixel in bad_pixels) or (tuple(next_pixel) in bad_pixels):
        # print("next_pixel in recursion", next_pixel)
        return right_consec_search(next_pixel, bad_pixels) + 1
    else:
        return 0


def consec_badpixels(bad_pixels):
    """Check for consecutive badpixels in the same nod.

    Consecutive in in axis=1 for the same axis=0 value.

    parameters
    ----------
    bad_pixels: list of list of ints
        List of index locations of bad pixels [nod, pixel].

    returns
    -------
    is_consec: bool
        True if a consecutve bad pixel found."""

    for pixel in bad_pixels:
        left_pix = [pixel[0], pixel[1] - 1]
        right_pix = [pixel[0], pixel[1] + 1]
        if (left_pix in bad_pixels) or (right_pix in bad_pixels):
            return True

    return False


def warn_consec_badpixels(bad_pixels, stop=True):
    """Raise erro on consecutive badpixels in the same nod.

    parameters
    ----------
    bad_pixels: list of list of ints
        List of index locations of bad pixels [nod, pixel]."""

    if consec_badpixels(bad_pixels):
        if stop:
            raise ValueError("Consective bad pixels were found. Need to deal with these.")
        else:
            logging.warning("Consective bad pixels in a nod were found.")
    return None
