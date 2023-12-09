#This script contains the functions used to transform spectra prior to computing similarity scores

import scipy.stats
import numpy as np

def wf_transform(wf_int, wf_mz, spec):
    #This function performs a weight factor transformation to a spectrum
    
    #input:
    #wf_int: float
    #wf_mz: float
    #spec: nx2 np array with first column being mass/charge and second column being intensity

    #output:
    #weight factor transformed spectrum

    spec[:,1] = np.power(spec[:,0], wf_mz) * np.power(spec[:,1], wf_int)
    return(spec)


def transform_int(intensity, thresh):
    #This transformation was presented by: 
    #Li, Y.; Kind, T.; Folz, J.; Vaniya, A.; Mehta, S. S.; Fiehn, O.
    #Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification. 
    #Nature Methods 2021, 18, 1524–1531
    
    #input:
    #intensity: 1d np array
    #thresh: nonnegative float
    
    #output:
    #1d np array of transformed intensities

    S = scipy.stats.entropy(intensity)
    if S < thresh:
        w = (1 + S) / (1 + thresh) 
        intensity = np.power(intensity, w)
    return intensity 


def normalize(intensity):
    #Normalizes a given vector to sum to 1 so that it represents a probability distribution
    intensity /= np.sum(intensity)
    return(intensity)


def clean_spectrum(spectrum, noise_removal, da):
    #This function was presented by: 
    #Li, Y.; Kind, T.; Folz, J.; Vaniya, A.; Mehta, S. S.; Fiehn, O.
    #Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification. 
    #Nature Methods 2021, 18, 1524–1531

    #This function: 
    #1) centroids peaks by merging peaks within a given window-size (i.e. da parameter)
    #2) removes peaks that have intensity lower than max(intensity)*noise_removal
    #3) normalizes the centroided and noise_removaled spectrum

    #input:
    #spectrum: nx2 np array with first column being mass/charge and second column being intensity
    #noise_removal and da parameters described above

    #output:
    #centroided, noise-removed, and normalized spectrum

    #Centroid peaks
    spectrum = spectrum[np.argsort(spectrum[:, 0])]
    spectrum = centroid_spec(spectrum, da=da)

    #Remove noise ions
    if noise_removal is not None and spectrum.shape[0] > 0:
        max_intensity = np.max(spectrum[:, 1])
        spectrum = spectrum[spectrum[:, 1] >= max_intensity * noise_removal]

    #Normalize the spectrum
    spectrum = normalize(spectrum)
    return spectrum


def centroid_spec(spec, da):
    #This function was presented by: 
    #Li, Y.; Kind, T.; Folz, J.; Vaniya, A.; Mehta, S. S.; Fiehn, O.
    #Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification. 
    #Nature Methods 2021, 18, 1524–1531

    #input:
    #spectrum: nx2 np array with first column being mass/charge and second column being intensity
    #da: window-size parameter

    #output:
    #centroided spectrum

    #Fast check is the spectrum needs centroiding
    mz_array = spec[:, 0]
    need_centroid = 0
    if mz_array.shape[0] > 1:
        mz_delta = mz_array[1:] - mz_array[:-1]
        if np.min(mz_delta) <= da:
            need_centroid = 1

    if need_centroid:
        intensity_order = np.argsort(-spec[:, 1])
        spec_new = []
        for i in intensity_order:
            mz_delta_allowed = da

            if spec[i, 1] > 0:
                #Find left board for current peak
                i_left = i - 1
                while i_left >= 0:
                    mz_delta_left = spec[i, 0] - spec[i_left, 0]
                    if mz_delta_left <= mz_delta_allowed:
                        i_left -= 1
                    else:
                        break
                i_left += 1

                #Find right board for current peak
                i_right = i + 1
                while i_right < spec.shape[0]:
                    mz_delta_right = spec[i_right, 0] - spec[i, 0]
                    if mz_delta_right <= mz_delta_allowed:
                        i_right += 1
                    else:
                        break

                #Merge those peaks
                intensity_sum = np.sum(spec[i_left:i_right, 1])
                intensity_weighted_sum = np.sum(spec[i_left:i_right, 0] * spec[i_left:i_right, 1])

                spec_new.append([intensity_weighted_sum / intensity_sum, intensity_sum])
                spec[i_left:i_right, 1] = 0

        spec_new = np.array(spec_new)
        spec_new = spec_new[np.argsort(spec_new[:, 0])]
        return spec_new
    else:
        return spec



def match_peaks_in_spectra(spec_a, spec_b, da):
    #This function was presented by: 
    #Li, Y.; Kind, T.; Folz, J.; Vaniya, A.; Mehta, S. S.; Fiehn, O.
    #Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification. 
    #Nature Methods 2021, 18, 1524–1531

    #This function matches two spectra to find common peaks in order
    #to obtain two lists of intensities of the same length

    #input:
    #spec_a: nx2 np array with first column being mass/charge and second column being intensity
    #spec_b: mx2 np array with first column being mass/charge and second column being intensity
    #da: window-size parameter

    #output:
    #kx3 np array with first column being mass/charge, second column being matched intensities of spec_a, and third column being matched intensities of spec_b

    a = 0
    b = 0

    spec_merged = []
    peak_b_int = 0.

    while a < spec_a.shape[0] and b < spec_b.shape[0]:
        mass_delta = spec_a[a, 0] - spec_b[b, 0]
        
        if mass_delta < -da:
            # Peak only existed in spec a.
            spec_merged.append([spec_a[a, 0], spec_a[a, 1], peak_b_int])
            peak_b_int = 0.
            a += 1
        elif mass_delta > da:
            # Peak only existed in spec b.
            spec_merged.append([spec_b[b, 0], 0., spec_b[b, 1]])
            b += 1
        else:
            # Peak existed in both spec.
            peak_b_int += spec_b[b, 1]
            b += 1

    if peak_b_int > 0.:
        spec_merged.append([spec_a[a, 0], spec_a[a, 1], peak_b_int])
        peak_b_int = 0.
        a += 1

    if b < spec_b.shape[0]:
        spec_merged += [[x[0], 0., x[1]] for x in spec_b[b:]]

    if a < spec_a.shape[0]:
        spec_merged += [[x[0], x[1], 0.] for x in spec_a[a:]]

    if spec_merged:
        spec_merged = np.array(spec_merged, dtype=np.float64)
    else:
        spec_merged = np.array([[0., 0., 0.]], dtype=np.float64)
    return spec_merged




