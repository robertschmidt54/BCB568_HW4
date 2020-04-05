# Imports
import numpy as np
from Bio import SeqIO
import pandas as pd

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

def alpha_func(initial_n, n, eta0, eta1, pi, gamma, reads, qualities, tau, h):

    alpha_vec = np.zeros(len(reads))
    nuc = ['A', 'C', 'G', 'T']

    list_of_hstrings = []
    sample_space = generate_hstring(list_of_hstrings, nuc, "", len(nuc), h)




    for i in range(0, len(reads)):

        q = qualities[i]
        r = reads[i]
        alpha = 1
        if r[0:h] == initial_n:
            # Is this Pi(R,N) or Pir(N,R)???

            for l in range(0, h):
                alpha *= eta0[q[l]] * pi[code[initial_n[l]]][code[r[l]]] * gamma[code[initial_n[l]]]
        else:
            for l in range(0, h):
                alpha *= eta1[q[l]] * pi[code[initial_n[l]]][code[r[l]]] * gamma[code[initial_n[l]]]

        for j in range(1,len(reads[i])):
            sum_thing = 0

            if r[j] == n:

                for k in nuc:
                    sum_thing += tau[k][code[n]] * alpha

                alpha = eta0[q[0]] * pi[code[n]][code[r[0]]] * sum_thing

            else:

                for k in nuc:
                    sum_thing += tau[k][code[n]] * alpha

                    alpha = eta1[q[0]] * pi[code[n]][code[r[0]]] * sum_thing
        alpha_vec[i] = alpha

    return alpha_vec







##### To Do: #####
# Eta0 and Eta1: Rob DONE
# Pi: Rob DONE
# Foward(alpha) : Kelby & Parnal
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
h = 3

# Import Fastq
reads = []
qualities = []
# Initial state matrix
tau = intialize_Tau(nucs, h)
gamma = np.random.dirichlet(np.ones(pow(4,h)))
eta0 = np.random.dirichlet(np.ones(40))
eta1 = np.random.dirichlet(np.ones(40))
# print(eta1)

pi = np.zeros((4,4))
for i in range(0,4):
    pi[i] = np.random.dirichlet(np.ones(4))


for record in SeqIO.parse("test.fastq", "fastq"):

    reads.append(str(record.seq))
    qualities.append(record.letter_annotations["phred_quality"])

# Convert to Numpy Arrays
reads = np.asarray(reads)
qualities = np.asarray(qualities)

print(len(reads))
print(alpha_func('TTA','G', eta0, eta1, pi, gamma, reads, qualities, tau, h))



# Export Data