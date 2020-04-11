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


def eijzw_num(alpha, beta, pi, eta0, eta1, tau, r, q, sample_space):


    # alpha is the j-1 row of the alpha matrix for 1 read
    # beta is the jth(last row = beta[1]) row of the beta matrix
    nucs = ['A','C','G','T']

    e_num = np.zeros((pow(4,h),4))
    # zloop
    for z in sample_space:
        # wloop
        temp_e = np.ones(4)
        for w in nucs:
            # Equation to calculate num of eijzw: alpha(z) * beta(w) * eta0(r=w) or eta1(r != w) * tau zw * pi(r,w)
            if r == w:
                temp_e[nucs.index(w)] = alpha[sample_space.index(z)] * beta[nucs.index(w)] * eta0[q-1] * tau[w][z] * pi[w][r]
            else:
                temp_e[nucs.index(w)] = alpha[sample_space.index(z)] * beta[nucs.index(w)] * eta1[q- 1] * tau[w][z] * pi[w][r]

        e_num[sample_space.index(z)] = temp_e

    sum_eijzw = np.sum(e_num)
    e = e_num/sum_eijzw
    return e




def e_func(eta0, eta1, pi, gamma, read, quality, tau, h):

    beta_out = np.ones((2, pow(4, h)))

    nuc = ['A', 'C', 'G', 'T']
    list_of_hstrings = []
    sample_space = generate_hstring(list_of_hstrings, nuc, "", len(nuc), h)
    alpha = alpha_func(eta0, eta1, pi, gamma, read, quality, tau, h)
    E_ijn = dict.fromkeys(sample_space, [])

    E_ijzw = np.zeros((pow(4, h), 4))
    e = np.zeros((len(read), pow(4,h)))
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
        e[j] = np.multiply(alpha[j], beta_out[1])


        # Calculate the numerator of E_ijzw for 1 read
        a = alpha[j-1]
        b = beta_out[1]
        E_ijzw += eijzw_num(a, b, pi, eta0, eta1, tau, read[j], quality[j], sample_space)
        beta_out[1] = beta_out[0]

    for z in sample_space:
        E_ijn[z] = e[:, sample_space.index(z)]

    eijn_sums = []
    for k in range(len(read)):
        position_sum = 0
        for keys in E_ijn.keys():
            position_sum += E_ijn[keys][k]

        eijn_sums.append(position_sum)


    # Dividing by the sum over all nucleotides within the same poisition same read
    for z in sample_space:
        E_ijn[z] = [b / m for b, m in zip(E_ijn[z], eijn_sums)]

    return E_ijn, E_ijzw


def Update_Gamma(ei1n, y, h):

    return np.sum(ei1n, axis=1)/y



def Update_Pi():
    pass

def Update_Eta():
    pass

def Update_Tau():
    pass




##### To Do: #####
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

# Import Fastq
reads = []
qualities = []

for record in SeqIO.parse("test.fastq", "fastq"):

    reads.append(str(record.seq))
    qualities.append(record.letter_annotations["phred_quality"])

# Convert to Numpy Arrays
reads = np.asarray(reads)
qualities = np.asarray(qualities)


# Da Real stuff
# Order of markov chain
h = 1

# Initial state matrix
tau = intialize_Tau(nucs, h)
gamma = np.random.dirichlet(np.ones(pow(4,h)))
eta0 = np.random.dirichlet(np.ones(40))
eta1 = np.random.dirichlet(np.ones(40))
pi = np.zeros((4,4))
list_of_hstrings = []
sample_space = generate_hstring(list_of_hstrings, nucs, "", len(nucs), h)
for i in range(0,4):
    pi[i] = np.random.dirichlet(np.ones(4))
pi = pd.DataFrame(pi, index = ['A', 'C', 'G', 'T'], columns = ['A', 'C', 'G', 'T'])

# # Test Crap
# eta0 = np.zeros(40)
# eta0[36:40] = [0.1, 0.1, 0.3, 0.4]
#
# eta1 = eta0 #np.random.dirichlet(np.ones(40))
# pi = np.array([[0.5, 0.25, 0.125, 0.125],[0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25]])
# pi = pd.DataFrame(pi, index= ('A', 'C', 'G', 'T'), columns= ('A', 'C', 'G', 'T'))
#
# tau = pd.DataFrame(np.array([[0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25]]), index=['A', 'C', 'G', 'T'], columns=['A', 'C', 'G', 'T'])
# gamma = np.array([0.9,0.03,0.04,0.03])
#
#
E_ijzw = np.zeros((pow(4, h),4))
E_ijn = np.ones((pow(4,h), len(reads),len(reads[0])))

for i in range(0,len(reads)):
    temp_eijn = e_func(eta0, eta1, pi, gamma, reads[i], qualities[i], tau, h)[0]

    E_ijzw += e_func(eta0, eta1, pi, gamma, reads[i], qualities[i], tau, h)[1]
    # print(E_ijn)
    for z in sample_space:
        E_ijn[sample_space.index(z)][i] = temp_eijn[z]

print("Eijn", E_ijn)
print(E_ijn[:,:,0])
gamma_new = Update_Gamma(E_ijn[:,:,0], len(reads), h)
# print(E_ijzw)







# Export Data