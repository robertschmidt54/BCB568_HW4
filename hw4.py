# Imports
import numpy as np
from Bio import SeqIO
import pandas as pd
import math

# Functions

def intialize_Tau(nucs, h):
    """

    :param nucs:
    :param h:
    :return:
    """
    list_of_hstrings = []
    generate_hstring(list_of_hstrings, nucs, "", len(nucs), h)

    matrix = np.random.rand(pow(4,h), 4)
    test = matrix / matrix.sum(axis=1)[:, None]

    Tau = pd.DataFrame(test, index= list_of_hstrings, columns= ('A', 'C', 'G', 'T'))

    return Tau


def generate_hstring(h_list, set_char, prefix, set_len, str_len):
    """

    :return:
    """

    if(str_len == 0):
        # Append to list
        h_list.append(prefix)
        return

    for i in range(set_len):
        # Next character of input added
        newPrefix = prefix + set_char[i]

        generate_hstring(h_list, set_char, newPrefix, set_len, str_len - 1)

    return h_list


def alpha_func(eta0, eta1, pi, gamma, read, quality, tau, h):

    alpha_vec = np.ones((len(read), pow(4,h)))
    nuc = ['A', 'C', 'G', 'T']
    list_of_hstrings = []
    sample_space = generate_hstring(list_of_hstrings, nuc, "", len(nuc), h)

    for s in sample_space:
        if read[0:h] == s:
            # Assume the first H-mer is error free
            for l in range(0, h):
                # Loop through first h-mer
                alpha_vec[0][sample_space.index(s)] *= eta0[quality[l]-1]*pi[s[l]][read[l]]*gamma[code[s[l]]]
        else:
            alpha_vec[0][sample_space.index(s)] = 0

    for j in range(h + 1, len(read)+1):
        for s in sample_space:
            alpha_vec[j-h][sample_space.index(s)] = 0

            for t in sample_space:
                alpha_vec[j-h][sample_space.index(s)] += alpha_vec[j-h-1][sample_space.index(t)]*tau[s[h-1]][t]

            if read[j-h] == s[h-1]:
                alpha_vec[j-h][sample_space.index(s)] *= eta0[quality[j-h]-1]*pi[s[h-1]][read[j-h]]
            else:
                alpha_vec[j-h][sample_space.index(s)] *= eta1[quality[j-h]-1] * pi[s[h-1]][read[j-h]]

    return alpha_vec


def e_func(eta0, eta1, pi, gamma, read, quality, tau, h):

    beta_out = np.ones((2, pow(4, h)))

    nuc = ['A', 'C', 'G', 'T']
    list_of_hstrings = []
    sample_space = generate_hstring(list_of_hstrings, nuc, "", len(nuc), h)
    alpha = alpha_func(eta0, eta1, pi, gamma, read, quality, tau, h)
    E_ijn = dict.fromkeys(['A', 'C', 'G', 'T'])
    E_ijn['A']= []
    E_ijn['C'] =[]
    E_ijn['G']= []
    E_ijn['T'] =[]

    E_ijzw = dict.fromkeys(['A', 'C', 'G', 'T'])
    E_ijn['A'] = []
    E_ijn['C'] = []
    E_ijn['G'] = []
    E_ijn['T'] = []

    for j in range(len(read) - 1, -1, -1):
        # Calculate beta

        for s in sample_space:
            beta_out[0][sample_space.index(s)] = 0
            for t in sample_space:

                # if read[j+1] == s[h-1]:
                if read[j] == s[h - 1]:
                    # beta_out[0][sample_space.index(s)] += beta_out[1][sample_space.index(t)] * eta0[quality[j+1]-1] * pi[t[h-1]][read[j+1]] * tau[t[h-1]][s]
                    beta_out[0][sample_space.index(s)] += beta_out[1][sample_space.index(t)] * eta0[
                        quality[j] - 1] * pi[t[h - 1]][read[j]] * tau[t[h - 1]][s]
                else:
                    # beta_out[0][sample_space.index(s)] += beta_out[1][sample_space.index(t)] * eta1[quality[j+1] - 1] * pi[t[h-1]][read[j+1]] * tau[t[h-1]][s]
                    beta_out[0][sample_space.index(s)] += beta_out[1][sample_space.index(t)] * eta1[
                        quality[j] - 1] * pi[t[h - 1]][read[j]] * tau[t[h - 1]][s]

        # Calculate the numerator of E_ijn for 1 read
        e = np.multiply(alpha[j], beta_out[1])
        E_ijn['A'].insert(0,e[0])
        E_ijn['C'].insert(0,e[1])
        E_ijn['G'].insert(0,e[2])
        E_ijn['T'].insert(0,e[3])

        # Calculate the numerator of E_ijzw for 1 read

        e_zw = np.multiply(alpha[j-1], beta_out[1])

        beta_out[1] = beta_out[0]

    eijn_sums = []
    for k in range(len(read)):
        position_sum = 0
        for keys in E_ijn.keys():
            position_sum += E_ijn[keys][k]

        eijn_sums.append(position_sum)

    # Dividng by the sum over all nucleotides within the same poisition.same read
    E_ijn['A'] = [b / m for b, m in zip(E_ijn['A'], eijn_sums)]
    E_ijn['C'] = [b / m for b, m in zip(E_ijn['C'], eijn_sums)]
    E_ijn['G'] = [b / m for b, m in zip(E_ijn['G'], eijn_sums)]
    E_ijn['T'] = [b / m for b, m in zip(E_ijn['T'], eijn_sums)]



    print(E_ijn)


    return E_ijn



    # Return:
    # ei1n,
    # sum over all j of eijn,
    # sum over all j and all nucs of eijk,
    # sum over all j






##### To Do: #####
# Estep:
    # eij:
    # eijz:
# Mstep:
    # Update Eta:
    # Update Pi:
    # Update Gamma
    # Update Tau
# Convergence function:
# Main Program:
    # Intialize Parameters
    # Run the Estep functions
    # Run Update functions

#to correct for underflow:
# calculate e numerators in log space
# Take the max of that numerator
# Convert back to normal with exp()
# Subtract the max value
# https://stackoverflow.com/questions/42599498/numercially-stable-softmax



code={"A":0, "C":1, "G":2, "T":3 }

nucs = ['A', 'C', 'G', 'T']
# Order of markov chain
h = 2

# Import Fastq
reads = []
qualities = []
# Initial state matrix
tau = intialize_Tau(nucs, h)
gamma = np.random.dirichlet(np.ones(pow(4,h)))
eta0 = np.random.dirichlet(np.ones(40))
eta1 = np.random.dirichlet(np.ones(40))


pi = np.zeros((4,4))
for i in range(0,4):
    pi[i] = np.random.dirichlet(np.ones(4))

pi = pd.DataFrame(pi, index = ['A', 'C', 'G', 'T'], columns = ['A', 'C', 'G', 'T'])
for record in SeqIO.parse("robs_test.fastq", "fastq"):

    reads.append(str(record.seq))
    qualities.append(record.letter_annotations["phred_quality"])

# Convert to Numpy Arrays
reads = np.asarray(reads)
qualities = np.asarray(qualities)

# h = 1
# eta0 = np.zeros(40)
# eta0[36:40] = [0.1, 0.1, 0.3, 0.4]
#
# eta1 = eta0 #np.random.dirichlet(np.ones(40))
# pi = np.array([[0.5, 0.25, 0.125, 0.125],[0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25]])
# pi = pd.DataFrame(pi, index= ('A', 'C', 'G', 'T'), columns= ('A', 'C', 'G', 'T'))
#
# tau = pd.DataFrame(np.array([[0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25]]), index=['A', 'C', 'G', 'T'], columns=['A', 'C', 'G', 'T'])
# gamma = np.array([0.9,0.03,0.04,0.03])
# alpha = alpha_func(eta0, eta1, pi, gamma, reads[0], qualities[0], tau, h)
# # print(alpha)
e_func(eta0, eta1, pi, gamma, reads[0], qualities[0], tau, h)
# print(alpha)


# Export Data