
def s_profile(wave, wav_0, depth, width, k):
    """Generate an S-profile from appendix of ferluga et al. """

    s_amp = 2 * depth * np.exp(-np.pi * depth**2 * ((wave - wav_0)**2 + (k/2)**2) / width**2)
    s_inner = (np.pi * depth**2 * (wave - wav_0) * k) / width**2
    s_lambda = s_amp * np.sinh(s_inner)
    return s_lambda


# I1, W1 and wav_1 are peak parameters.
