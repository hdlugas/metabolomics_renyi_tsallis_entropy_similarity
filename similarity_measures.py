##### Similarity Score Functions #####

import scipy.stats
import numpy as np


def S_cos(X, Y):
    #Cosine Similarity Measure
    return np.dot(X,Y) / (np.sqrt(sum(np.power(X,2))) * np.sqrt(sum(np.power(Y,2))))

def S_shannon(ints_a, ints_b):
    #Shannon Entropy Similarity Measure
    #This similarity function was presented by: 
    #Li, Y.; Kind, T.; Folz, J.; Vaniya, A.; Mehta, S. S.; Fiehn, O.
    #Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification. 
    #Nature Methods 2021, 18, 1524–1531
    #Since scipy.stats.entropy normalizes the input vector to sum to 1, vec1 and vec1 need not be normalized
    ent_a = scipy.stats.entropy(ints_a)
    ent_b = scipy.stats.entropy(ints_b)
    ent_ab = scipy.stats.entropy(ints_a + ints_b)
    return(1 - (2 * ent_ab - ent_a - ent_b)/np.log(4))

def ent_renyi(ints, q):
    #Computes the Renyi entropy of a probability distribution for a given positive entropy dimension q
    return np.log(sum(np.power(ints,q))) / (1-q)

def ent_tsallis(ints, q):
    #Computes the Tsallis entropy of a probability distribution for a given positive entropy dimension q
    return (sum(np.power(ints,q))-1) / (1-q)

def S_renyi(ints_a, ints_b, q):
    #Renyi Entropy Similarity Measure
    #This is a novel similarity measure which generalizes the Shannon Entropy Similarity Measure
    #Note that the Renyi Similarity Measure approaches the Shannon Entropy Similiarity Measure as q approaches 1
    #Note that ints_a and ints_b must be normalized to sum to 1
    ent_a = ent_renyi(ints_a, q)
    ent_b = ent_renyi(ints_b, q)
    ent_merg = ent_renyi(ints_a/2 + ints_b/2, q)
    N = (1/(1-q)) * (2*np.log(np.sum(np.power(ints_a/2,q))+np.sum(np.power(ints_b/2,q))) - np.log(np.sum(np.power(ints_a,q))) - np.log(np.sum(np.power(ints_b,q))))
    return 1 - (2 * ent_merg - ent_a - ent_b) / N

def S_tsallis(ints_a, ints_b, q):
    #Tsallis Entropy Similarity Measure
    #This is a novel similarity measure which generalizes the Shannon Entropy Similarity Measure
    #Note that the Tsallis Similarity Measure approaches the Shannon Entropy Similiarity Measure as q approaches 1
    #Note that ints_a and ints_b must be normalized to sum to 1
    ent_a = ent_tsallis(ints_a, q)
    ent_b = ent_tsallis(ints_b, q)
    ent_merg = ent_tsallis(ints_a/2 + ints_b/2, q)
    N = np.sum(2*np.power(ints_a/2,q)+2*np.power(ints_b/2,q)-np.power(ints_a,q)-np.power(ints_b,q)) / (1-q)
    return 1 - (2 * ent_merg - ent_a - ent_b) / N




