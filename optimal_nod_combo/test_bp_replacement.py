"""Tests for bad pixel replacement."""
import pytest
import numpy as np
import bp_replacement as pb


def test_sigma_detect():
    """Create a dataset. Add a value that is outside 4 sigma.

    See if they are identified correctly.
    """
    data = np.random.normal(loc=1, scale=.5, size=(8, 50))

    bad_pixels = [(4, 1), (3, 10), (7, 25), (2, 40), (4, 48)]
    for bad_loc in bad_pixels:
        if bad_loc[1] < 2:
            slice_data = data[:, :6]
        elif bad_loc[1] > data.shape[1]:
            slice_data = data[:, -5:]
        else:
            slice_data = data[:, bad_loc[1] - 2: bad_loc[1] + 3]
        data[bad_loc[0], bad_loc[1]] += (8 * np.std(slice_data.ravel()) *
                                         (2 * np.random.randint(0, 2) - 1))   # add sign

    bad_pixel_record = pb.sigma_detect(data, plot=False)

    # Test all pixels recovered
    assert len(bad_pixel_record) == len(bad_pixels)
    print("here i am")
    # Test the index values match
    assert np.all([pix_record[0] == bad_pixel[0]
                   for pix_record, bad_pixel in zip(bad_pixel_record, bad_pixels)])
    assert np.all([pix_record[1] == bad_pixel[1]
                   for pix_record, bad_pixel in zip(bad_pixel_record, bad_pixels)])


@pytest.mark.parametrize("test_input,expected", [
    ([[1, 10], [1, 11]], True),
    ([[5, 1], [5, 2]], True),
    ([[1, 10], [2, 10], [3, 10]], False),
    ([[1, 7], [2, 8], [4, 9]], False)])
def test_consecutive_badpixels(test_input, expected):
    """Test detection of consecutives."""
    assert pb.consec_badpixels(test_input) == expected


def test_warn_consec_badpixels():
    """Test warn consutive throws and error."""
    with pytest.raises(ValueError):
        pb.warn_consec_badpixels([[1, 10], [1, 11]])

    assert pb.warn_consec_badpixels([[1, 7], [2, 7], [3, 8]]) is None


def test_left_consecutive_count():
    """Test consecutive counting works on left."""
    bad_pixels = [[1, 2], [1, 3], [1, 4], [1, 5]]
    pixel = [1, 5]
    assert pb.left_consec_search(pixel, bad_pixels) == 3


def test_right_consecutive_count():
    """Test consecutive counting works on right."""
    bad_pixels = [[1, 2], [0, 8], [1, 3], [0, 0], [1, 4], [1, 5], [1, 7]]
    pixel = [1, 4]
    assert pb.right_consec_search(pixel, bad_pixels) == 1


def test_consecutive_counts():
    """Test consecutive counting works on both sides, and sum is correct."""
    bad_pixels = [[1, 2], [0, 8], [1, 3], [0, 0], [1, 4], [1, 5], [1, 8],
                  [2, 9], [1, 7], [1, 11], [1, 12], [1, 13]]
    pixel = [1, 6]
    assert pb.right_consec_search(pixel, bad_pixels) == 2
    assert pb.left_consec_search(pixel, bad_pixels) == 4
    assert (pb.right_consec_search(pixel, bad_pixels) +
            pb.left_consec_search(pixel, bad_pixels)) == 6


def test_interp_badpixel():
    """Test of Linear interpolation across bad pixels.

    Single bad pixels at the start, middle and end.
    """
    nods = np.array([[-99, 5, -99, 3, 1], [1, -99, 4, 1, 2], [1, 3, 4, 5, -99]], dtype=np.float32)
    interp_nods = pb.interp_badpixels(nods, [[0, 0], [0, 2], [1, 1], [2, 4]])

    expected_array = np.array([[5, 5, 4, 3, 1], [1, 2.5, 4, 1, 2], [1, 3, 4, 5, 5]])

    assert np.allclose(interp_nods, expected_array)


def test_bad_pixel_interp_types():
    """Nod must be array, bad_pixels must be list.

    interp_badpixels does a core dump if given a bad_pixels array.
    """
    bad_pixels = np.array([[1, 2], [3, 4]])

    with pytest.raises(TypeError):
        pb.interp_badpixels(np.array([1, 2, 3]), bad_pixels)   # second value should be a list

    with pytest.raises(TypeError):
        pb.interp_badpixels([1, 2, 3], list(bad_pixels))  # first values should not be a list


# Add in a number of examples. The first 9 are systematic. The last two are more complicated.
@pytest.mark.parametrize("test_in,bad,expected", [
    ([[-99, 2, 3, 4]], [[0, 0]], [[2, 2, 3, 4]]),
    ([[1, -99, 3, 4]], [[0, 1]], [[1, 2, 3, 4]]),
    ([[1, 2, -99, 4]], [[0, 2]], [[1, 2, 3, 4]]),
    ([[1, 2, 3, -99]], [[0, 3]], [[1, 2, 3, 3]]),
    ([[-99, -99, 3, 4]], [[0, 0], [0, 1]], [[3, 3, 3, 4]]),
    ([[1, -99, -99, 4]], [[0, 1], [0, 2]], [[1, 2, 3, 4]]),
    ([[1, 2, -99, -99]], [[0, 2], [0, 3]], [[1, 2, 2, 2]]),
    ([[-99, -99, -99, 4]], [[0, 1], [0, 2], [0, 0]], [[4, 4, 4, 4]]),
    ([[1, -99, -99, -99]], [[0, 2], [0, 3], [0, 1]], [[1, 1, 1, 1]]),
    ([[8, -99, -99, 5]], [[0, 1], [0, 2]], [[8, 7, 6, 5]]),
    ([[-99, -99, -99, 1.5, 3, 1, 3, 4, -99, 0], [3, 5, 1, 3, 4, -99, -99, -99, -99, 0]],
     [[0, 0], [0, 1], [0, 2], [0, 8], [1, 5], [1, 6], [1, 7], [1, 8]],
     [[1.5, 1.5, 1.5, 1.5, 3, 1, 3, 4, 2, 0], [3, 5, 1, 3, 4, 3.2, 2.4, 1.6, 0.8, 0]]),
    ([[1, 3, 4, -99, -99, 1, 3, 4, -99, 2], [3, 5, 1, 3, 4, -99, -99, -99, -99, 2]],
     [[0, 3], [0, 4], [0, 8], [1, 5], [1, 6], [1, 7], [1, 8]],
     [[1, 3, 4, 3, 2, 1, 3, 4, 3, 2], [3, 5, 1, 3, 4, 3.6, 3.2, 2.8, 2.4, 2]])
])
def test_multiple_bp_interpolation(test_in, bad, expected):
    """Linear interpolation with multiple bad pixels.

    Parametrized to catch all the edge cases. Using -99 for clarity of bad pixels.
    """
    nods = np.array(test_in, dtype=np.float32)
    expected_array = np.array(expected, dtype=np.float32)

    interp_nods = pb.interp_badpixels(nods, bad)
    print("input", nods)
    print("expected", expected_array)
    print("interped", interp_nods)

    assert np.allclose(interp_nods, expected_array)


@pytest.mark.parametrize("test_in", [
    ([[-99], [1], [2]]),
    ([[1, 2], [3, 4]])])
def test_fail_on_small_array(test_in):
    """Test small arrays are failed."""
    nods = np.array(test_in, dtype=np.float32)
    bad = [[0, 0], [1, 0]]

    with pytest.raises(ValueError):
        interp_nods = pb.interp_badpixels(nods, bad)


def test_allbadpixels():
    """A completly bad row is unchanged."""
    nods = np.array([[-99, -99, -99], [-99, -99, 1]], dtype=np.float32)
    expected_array = np.array([[-99, -99, -99], [1, 1, 1]], dtype=np.float32)

    bad = [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1]]

    interp_nods = pb.interp_badpixels(nods, bad)
    print("input", nods)
    print("expected", expected_array)
    print("interped", interp_nods)

    assert np.allclose(interp_nods, expected_array)
