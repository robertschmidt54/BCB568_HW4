
import numpy as np
from Bio import SeqIO
import pandas as pd
import csv

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



def Verterbi(pie, gamma, tau, eta0, eta1, reads, qualities, h, sample_space, read_names):


    true_seq_vec = []
    with open('corrected_error_final_data_h2.csv', mode = 'w') as fh:
        error_writer = csv.writer(fh, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        error_writer.writerow(["Name of Read", "Position of Error", "Read Base", "Quality Score of Base", "Correct Base"  ])
        for i in range(len(reads)):
            # Forward Process:
            read = reads[i]
            quality = qualities[i]
            Zeta = np.zeros((2, pow(4,h)))
            B = np.empty((len(reads[0]), pow(4, h)), dtype= object)
            true_seq = ""

            z = read[0:h]
            Zeta[0][sample_space.index(z)] = 1
            # Assume the first H-mer is error free
            for l in range(0, h):
                # Loop through first h-mer
                Zeta[0][sample_space.index(z)] *= eta0[quality[0]-1] * pie[z[l]][read[l]]
            Zeta[0][sample_space.index(z)] *= gamma[sample_space.index(z)]

            for j in range(h, len(read)):
                    # Calculate new Zetas = eta0 or eta1 * pie * max(Zeta1(z') * tau(z', z)

                    for z in sample_space:

                        if read[j] == z[h-1]:
                            Zeta[1][sample_space.index(z)] = eta0[quality[j] - 1] * pie[z[-1]][read[j]] * np.max((Zeta[0] * tau[z]).to_numpy())

                        else:
                            Zeta[1][sample_space.index(z)] = eta1[quality[j] - 1] * pie[z[-1]][read[j]] * np.max((Zeta[0] * tau[z]).to_numpy())

                        # Calculate a B matrix: Bj(z') = argmax(zeta_j-1(z') * tau(z', z) for j > 1
                        B[j][sample_space.index(z)] = sample_space[np.argmax((Zeta[0] * tau[z]).to_numpy())]

                    Zeta[0] = Zeta[1]

            for j in range(len(read)-1, h-2, -1):
                # Backward Process:

                if j == len(read)-1:
                    # Calculate last position in true sequence. argmax(gamma_l(z))
                    Z_hat = sample_space[np.argmax(Zeta[1])]
                    true_seq += Z_hat[-1]

                elif j == h-1:
                    new_Zhat = str(B[j + 1][sample_space.index(Z_hat)])
                    rev_Zhat = new_Zhat[::-1]
                    true_seq += rev_Zhat

                else:
                    new_Zhat = str(B[j+1][sample_space.index(Z_hat)])
                    new_nuc = new_Zhat[h-1]
                    true_seq += new_nuc
                    Z_hat = new_Zhat


            # Reverse the sequence
            true_seq_vec.append(true_seq[::-1])

        # Export Results
        for i in range(len(true_seq_vec)):
            for j in range(len(true_seq_vec[0])):

                if reads[i][j] == true_seq_vec[i][j]:
                    pass
                else:
                    error_writer.writerow([f'{read_names[i]}', f'{j}', f'{reads[i][j]}', f'{qualities[i][j]}',f'{true_seq_vec[i][j]}' ])

    fh.close()

    return


nucs = ['A','C','G','T']
h = 2
reads = []
qualities = []
read_names = []
for record in SeqIO.parse("errored_reads.fastq", "fastq"):
    reads.append(str(record.seq))
    read_names.append(str(record.id))
    qualities.append(record.letter_annotations["phred_quality"])
print("Data read in.")
# Convert to Numpy Arrays
reads = np.asarray(reads)
qualities = np.asarray(qualities)
list_of_hstrings = []
sample_space = generate_hstring(list_of_hstrings, nucs, "", len(nucs), h)

# When H=2
eta0 = np.array([0.00000000e+00, 7.48732417e-02, 0.00000000e+00, 7.24526631e-05,
                 2.94768134e-04, 3.87465277e-04, 4.96076520e-04, 4.62494727e-04,
                 7.74696116e-04, 6.94777618e-04, 4.61513277e-04, 5.59140658e-04,
                 7.42689260e-04, 2.05002015e-03, 2.43988494e-03, 1.40820910e-03,
                 1.28981445e-03, 1.72419227e-03, 1.93139254e-03, 2.33904544e-03,
                 3.31726665e-03, 2.95317448e-03, 4.46119144e-03, 4.37658945e-03,
                 6.33330969e-03, 6.32423606e-03, 8.34829415e-03, 1.14804003e-02,
                 2.02238679e-02, 2.19094548e-02, 4.09402674e-02, 5.37754488e-02,
                 1.24592137e-01, 5.38722638e-01, 5.92398484e-02, 0.00000000e+00,
                 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00])
eta1 = np.array([0.00000000e+00, 7.40514592e-02, 0.00000000e+00, 7.17765823e-05,
                 2.91102365e-04, 3.84335963e-04, 4.89219043e-04, 4.60882509e-04,
                 7.63140502e-04, 6.93790948e-04, 4.66809200e-04, 5.58440054e-04,
                 7.45507254e-04, 2.04766628e-03, 2.42561456e-03, 1.40394836e-03,
                 1.30785335e-03, 1.73862278e-03, 1.92921646e-03, 2.34063555e-03,
                 3.33481794e-03, 2.92606083e-03, 4.47139993e-03, 4.39796528e-03,
                 6.31894808e-03, 6.32750542e-03, 8.37239629e-03, 1.14982405e-02,
                 2.02748975e-02, 2.19363481e-02, 4.11298886e-02, 5.38083177e-02,
                 1.25241939e-01, 5.38312387e-01, 5.94788673e-02, 0.00000000e+00,
                 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00])

pi = np.array([[0.307758,  0.307826,  0.307808,  0.307868],
               [0.181251,  0.181192,  0.181232,  0.181218],
               [0.172631,  0.172585,  0.172591,  0.172574],
               [0.338360,  0.338397,  0.338369,  0.338339]])

pi = pd.DataFrame(pi, index = ['A', 'C', 'G', 'T'], columns = ['A', 'C', 'G', 'T'])

tau = np.array([[0.303195,  0.214504,  0.365723,  0.116578],
                [0.355689,  0.111196,  0.283540,  0.249574],
                [0.025283,  0.053001,  0.537225,  0.384491],
                [0.124580,  0.447810,  0.338510,  0.089100],
                [0.242212,  0.317849,  0.209377,  0.230562],
                [0.382567,  0.030090,  0.307903,  0.279440],
                [0.212504,  0.342151,  0.290727,  0.154618],
                [0.208341,  0.317893,  0.231864,  0.241902],
                [0.105379,  0.215727,  0.302115,  0.376779],
                [0.462467,  0.237296,  0.284047,  0.016190],
                [0.084682,  0.394116,  0.260716,  0.260486],
                [0.169698,  0.539953,  0.251723,  0.038626],
                [0.329776,  0.245832,  0.126980,  0.297412],
                [0.378216,  0.063844,  0.083876,  0.474064],
                [0.013416,  0.802889,  0.117358,  0.066337],
                [0.511813,  0.210880,  0.044435,  0.232872]])


new_tau = np.zeros((pow(4,h), pow(4,h)))
for i in sample_space:
    # Append A value
    next_state_A = i[1:h] + 'A'
    new_tau[sample_space.index(i)][sample_space.index(next_state_A)] = tau[sample_space.index(i)][0]

    next_state_C = i[1:h] + 'C'
    new_tau[sample_space.index(i)][sample_space.index(next_state_C)] = tau[sample_space.index(i)][1]

    next_state_G = i[1:h] + 'G'
    new_tau[sample_space.index(i)][sample_space.index(next_state_G)] = tau[sample_space.index(i)][2]

    next_state_T = i[1:h] + 'T'
    new_tau[sample_space.index(i)][sample_space.index(next_state_T)] = tau[sample_space.index(i)][3]


tau = pd.DataFrame(new_tau, index = sample_space, columns = sample_space)



gamma = np.array([0.06192457, 0.06927721, 0.06327057, 0.05552765, 0.06192457, 0.06927721,
                  0.06327057, 0.05552765, 0.06192457, 0.06927721, 0.06327057, 0.05552765,
                  0.06192457, 0.06927721, 0.06327057, 0.05552765])


Verterbi(pi, gamma, tau, eta0, eta1, reads, qualities, h, sample_space, read_names)



# # # WHEN H= 1
# eta0 = np.array([0.00000000e+00, 7.23827059e-02, 0.00000000e+00, 8.82718212e-05,
#                  3.09424990e-04, 4.21264276e-04, 5.04430509e-04, 4.65567974e-04,
#                  8.20572503e-04, 7.08005357e-04, 4.75448005e-04, 5.92311432e-04,
#                  7.68808474e-04, 2.07435052e-03, 2.46725869e-03, 1.44957936e-03,
#                  1.38155080e-03, 1.78990541e-03, 2.06269220e-03, 2.49898916e-03,
#                  3.47934889e-03, 3.08440625e-03, 4.57182739e-03, 4.59554365e-03,
#                  6.65917250e-03, 6.63957179e-03, 8.67933832e-03, 1.21120492e-02,
#                  2.13252629e-02, 2.21963374e-02, 4.21350498e-02, 5.51852124e-02,
#                  1.24508871e-01, 5.35307012e-01, 5.82598595e-02, 0.00000000e+00,
#                  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00])
# eta1 = np.array([0.00000000e+00, 7.48946385e-02, 0.00000000e+00, 6.63411547e-05,
#                  2.86025165e-04, 3.72699899e-04, 4.86260649e-04, 4.59806505e-04,
#                  7.47247612e-04, 6.89240223e-04, 4.62100998e-04, 5.47052355e-04,
#                  7.36585999e-04, 2.03928863e-03, 2.41603169e-03, 1.38969989e-03,
#                  1.27663016e-03, 1.71627721e-03, 1.88414868e-03, 2.28579329e-03,
#                  3.27945815e-03, 2.88067882e-03, 4.43358783e-03, 4.32314907e-03,
#                  6.20697646e-03, 6.21938342e-03, 8.25916796e-03, 1.12818142e-02,
#                  1.98977884e-02, 2.18383058e-02, 4.07226216e-02, 5.33251864e-02,
#                  1.25279296e-01, 5.39478460e-01, 5.98182575e-02, 0.00000000e+00,
#                  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00])
#
# pi = np.array([[0.307859,  0.307779,  0.307758,  0.307869],
#                [0.181122,  0.181255,  0.181361,  0.181236],
#                [0.172568,  0.172630,  0.172612,  0.172591],
#                [0.338451,  0.338336,  0.338269,  0.338304]])
#
#
# # pi = np.random.rand(4,4)
# # pi /= np.sum(pi, axis=0)
# pi = pd.DataFrame(pi, index = ['A', 'C', 'G', 'T'], columns = ['A', 'C', 'G', 'T'])
#
# #
# tau = np.array([[0.704595,  0.159773,  0.105305,  0.030328],
#                 [0.216565,  0.235261,  0.166677,  0.381497],
#                 [0.216266,  0.228128,  0.348351,  0.207254],
#                 [0.171331,  0.174885,  0.533431,  0.120352]])
#
# new_tau = np.zeros((pow(4,h), pow(4,h)))
# for i in sample_space:
#     # Append A value
#     next_state_A = i[1:h] + 'A'
#     new_tau[sample_space.index(i)][sample_space.index(next_state_A)] = tau[sample_space.index(i)][0]
#
#     next_state_C = i[1:h] + 'C'
#     new_tau[sample_space.index(i)][sample_space.index(next_state_C)] = tau[sample_space.index(i)][1]
#
#     next_state_G = i[1:h] + 'G'
#     new_tau[sample_space.index(i)][sample_space.index(next_state_G)] = tau[sample_space.index(i)][2]
#
#     next_state_T = i[1:h] + 'T'
#     new_tau[sample_space.index(i)][sample_space.index(next_state_T)] = tau[sample_space.index(i)][3]
# #
#
#
#
#
# tau = pd.DataFrame(new_tau, index = sample_space, columns = sample_space)
# gamma = np.array([0.41013573, 0.19297337, 0.24200037, 0.15489052])
#
#
# Verterbi(pi, gamma, tau, eta0, eta1, reads, qualities, h, sample_space, read_names)
#
#
#
#
#
#
#
