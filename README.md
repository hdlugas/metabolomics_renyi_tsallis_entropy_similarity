# Metabolomics: Compound Identification via Spectral Library Matching with Novel Entropy Similarity Measures

We have constructed two novel similarity measures based on Renyi and Tsallis Entropy to quantify the similarity between two probability distributions. Our original motivation for constructing these similarity measures is to improve compound identification via spectral library matching with respect to either liquid-chromatography mass spectrometry (LCMS) or gas-chromatography mass spectrometry (GCMS). We constructed these similarity measures by generalizing the Shannon Entropy Similarity Measure presented by Li et. al. in the publication referenced below.

Li, Y., Kind, T., Folz, J. et al. Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification. Nat Methods 18, 1524–1531 (2021). https://doi.org/10.1038/s41592-021-01331-z

We provide R and python implementations of the four similarity measures used in our analysis: Cosine, Shannon Entropy, Renyi Entropy, and Tsallis Entropy.
We also provide a python script with all functions used to process spectra prior to computing similarity scores. Note that some of these functions, along with the Shannon Entropy Similarity Measure, were originally presented by Li et. al., and we reference these functions in the scripts.

