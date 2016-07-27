import py.test


from IP_Convolution import  fast_wav_selector, unitary_Gauss, fast_convolve, IPconvolution

import numpy as np

def test_wav_selector():
    x = [1, 2, 3, 4,]
    y = [1, 2, 1, 1, 2]
    x_np = np.array(x)
    y_np = np.array(y)
    assert type(fast_wav_selector(x, y, 1, 3)) == list
    assert type(fast_wav_selector(x_np, y_np, 1, 3)) == list
    assert type(fast_wav_selector(x, y, 1, 3)[0]) == list
    assert type(fast_wav_selector(x, y, 1, 3)[1]) == list
    assert type(fast_wav_selector(x_np, y_np, 1, 3)[0]) == np.ndarray
    assert type(fast_wav_selector(x_np, y_np, 1, 3)[1]) == np.ndarray
    assert fast_wav_selector(x, y, 0, 3)[0] == [1,2]
    
#def test_unitary_Gauss():
#    FWHM = 1
#    normalization = 2 * np.sqrt(2 * np.log(2)) / (np.abs(FWHM)*np.sqrt(2*np.pi))
#    
#    assert unitary_Gauss(0, 0, FWHM) == normalization

    #assert type(unitary_Gauss([-1, 0, 1, 3], 0, FWHM)) == 'numpy.ndarray'


def test_fast_convolution():
    a = np.linspace(2130, 2170, 1024)
    b = np.linspace(2100, 2200, 1024)
    c = np.ones_like(b)
    R = 50000
    for a_val in a:
        #print(type(fast_convolve(a_val, R, b, c, 5)))
        #print(isinstance(fast_convolve(a_val, R, b, c, 5), np.float64))
        assert type(fast_convolve(a_val, R, b, c, 5)) == np.float64
        assert fast_convolve(a_val, R, b, c, 5) == 1     # Test a flat input of 1s gives a flat ouput of 1s
        assert fast_convolve(a_val, R, b, 0*c, 5) == 0     # Test a flat input of 1s gives a flat ouput of 1s




def test_IPconvolution():
	wave = [1,2,3,5,6,7,8,9,10,11] 
	flux = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
	chip_limits= [2,9]
	R = 100
	ans = IPconvolution(wave, flux, chip_limits, R, plot=False)
	assert ans == [1,1]












if __name__ == "__main__":

    test_fast_convolution()