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
h = 7

# Import Fastq
reads = []
qualities = []
# Initial state matrix
tau = intialize_Tau(nucs, h)
gamma = np.random.dirichlet(np.ones(pow(4,h)))
eta0 = np.random.dirichlet(np.ones(40))
eta1 = np.random.dirichlet(np.ones(40))
pi = np.random.dirichlet(np.ones((4,4)),4)
for record in SeqIO.parse("test.fastq", "fastq"):

    reads.append(str(record.seq))
    qualities.append(record.letter_annotations["phred_quality"])

# Convert to Numpy Arrays
reads = np.asarray(reads)
qualities = np.asarray(qualities)



# Export Data