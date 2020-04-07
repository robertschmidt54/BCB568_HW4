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

    alpha_vec = np.zeros((len(reads[0]), pow(4,h)))
    nuc = ['A', 'C', 'G', 'T']
    list_of_hstrings = []
    sample_space = generate_hstring(list_of_hstrings, nuc, "", len(nuc), h)
    position = 0
    j = 0

    # while j <= len(read):
    #     # print(j)
    #     if j > 0:
    #         # Correct for Underflow
    #
    #         max_value = np.max(alpha_vec[position-1])
    #         # alpha_vec[position - 1] = alpha_vec[position - 1] - max_value
    #         alpha_vec[position-1] = np.exp(alpha_vec[position-1])
    #         alpha_vec[position-1] = alpha_vec[position-1] - max_value
    #         alpha_vec[position-1] = alpha_vec[position-1]/np.sum(alpha_vec[position-1])
    #
    #         if j == len(reads[0]):
    #             break
        # print(j)
        # print(alpha_vec)

    for s in sample_space:
        alpha = 1

        if read[0:h] == s:
            # Assume the first H-mer is error free
        # if j == 0:
            for l in range(0, h):

                # Loop through first h-mer

                alpha_vec[0][sample_space.index(s)] += math.log(eta0[quality[l]-1] * pi[read[l]][s[l]] * gamma[code[s[l]]])

                # alpha *= eta0[quality[l] - 1] * pi[code[read[l]]][code[s[l]]] * gamma[code[s[l]]]
        else:
            alpha_vec[0][sample_space.index(s)] = 0

        for j in range(h + 1, len(read)):
            # for s in sample_space:
            for s in nuc:
                alpha_vec[j-h][sample_space.index(s)] = 0

                for t in sample_space:
                    alpha_vec[j - h][sample_space.index(s)] += math.exp(alpha_vec[j-h-1][sample_space.index(t)]) * tau[t][s]

                    # alpha[j - h][s] += alpha[j - h - 1][t] * tau[t][s]

            if read[j] == s[h-1]:

                alpha_vec[j - h][sample_space.index(s)] = math.log(alpha_vec[j - h][sample_space.index(s)]) + math.log(eta0[quality[j]-1] * pi[read[j]][s[h-1]])
                # n = s[h-1]
                # sum_thing = tau.sum(axis=0)[n]
                # # Something is wrong here.
                # sum_thing *= alpha_vec[position-1][sample_space.index(s)]
                # alpha += math.log(eta0[q[j]-1] * pi[code[r[j]]][code[n]] * sum_thing)

            else:

                alpha_vec[j - h][sample_space.index(s)] = math.log(alpha_vec[j - h][sample_space.index(s)]) + math.log(eta1[quality[j]-1] * pi[read[j]][s[h-1]])
                # n = s[h - 1]
                # sum_thing = tau.sum(axis=0)[n]
                # sum_thing *= alpha_vec[position - 1][sample_space.index(s)]
                # alpha += math.log(eta1[q[j]-1] * pi[code[r[j]]][code[n]] * sum_thing)
                # print(math.log(eta1[q[j]-1] * pi[code[r[j]]][code[n]] * sum_thing))

    return alpha_vec


def beta_func(eta0, eta1, pi, gamma, reads, qualities, tau, h):
    beta_out = np.ones((len(reads[0]), pow(4, h)))
    beta_l = 1

    for j in range(len(reads[0])-1, -1, -1):
        for i in len(reads):
            read = reads[i]
            quality = qualities[i]
            beta_ijn = 0
            for n in nucs:
                if read[j] == n:
                    beta_ijn += eta0[quality]*pi[1]






##### To Do: #####
# Eta0 and Eta1: Rob DONE
# Pi: Rob DONE
# Foward(alpha) : Kelby & Parnal DONE
# Backward(beta): Shatabdi & ROb
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
h = 7

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


for record in SeqIO.parse("robs_test.fastq", "fastq"):

    reads.append(str(record.seq))
    qualities.append(record.letter_annotations["phred_quality"])

# Convert to Numpy Arrays
reads = np.asarray(reads)
qualities = np.asarray(qualities)

h = 1
eta0 = np.zeros(40)
eta0[36:40] = [0.1, 0.1, 0.3, 0.4]

eta1 = np.random.dirichlet(np.ones(40))
pi = np.array([[0.5, 0.25, 0.125, 0.125],[0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25]])
pi = pd.DataFrame(pi, index= ('A', 'C', 'G', 'T'), columns= ('A', 'C', 'G', 'T'))

tau = intialize_Tau(nucs, h)
gamma = np.array([0.9,0.03,0.04,0.03])
alpha = alpha_func(eta0, eta1, pi, gamma, reads[0], qualities[0], tau, h)
print(alpha)


# Export Data