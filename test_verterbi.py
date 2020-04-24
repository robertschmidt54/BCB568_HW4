
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
    with open('corrected_error.csv', mode = 'w') as fh:
        error_writer = csv.writer(fh, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        error_writer.writerow(["Name of Read", "Position of Error", "Read Base", "Quality Score of Base", "Correct Base"  ])
        for i in range(len(reads)):

        # Forward Part:
            read = reads[i]
            # print("Read", read)
            quality = qualities[i]
            Zeta = np.zeros((2, pow(4,h)))
            # B = np.zeros((len(reads[0]),pow(4,h)))
            B = np.empty((len(reads[0]), pow(4, h)), dtype= object)
            true_seq = ""
            # for j in range(0, len(read)):
            # Forward Process
                # if j == 0:
            # for z in nucs:

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

            # print("B:", B)
            for j in range(len(read)-1, h-2, -1):
                # print("J", j)
                # print(read[j])
                # Backward Process:
                # THis works for h=1
                if j == len(read)-1:
                    # Calculate last position in true sequence. argmax(gamma_l(z))
                    # for z in sample_space:
                    # Z_hat = sample_space[np.argmax((Zeta[0] * tau[z[-1]]).to_numpy())]
                    Z_hat = sample_space[np.argmax(Zeta[1])]

                    # rev_Zhat = Z_hat[::-1]
                    true_seq += Z_hat[-1]

                elif j == h-1:
                    new_Zhat = str(B[j + 1][sample_space.index(Z_hat)])
                    rev_Zhat = new_Zhat[::-1]
                    true_seq += rev_Zhat
                    # true_seq += new_Zhat
                else:
                    # Zj = Bj+1(Zj+1)
                    new_Zhat = str(B[j+1][sample_space.index(Z_hat)])
                    new_nuc = new_Zhat[0]
                    true_seq += new_nuc
                    Z_hat = new_Zhat


            # Reverse the sequence
            # print("Read seq", read)
            # print("True seq", true_seq[::-1])
            true_seq_vec.append(true_seq[::-1])
        for i in range(len(true_seq_vec)):
            for j in range(len(true_seq_vec[0])):

                if reads[i][j] == true_seq_vec[i][j]:
                    pass
                else:
                    error_writer.writerow([f'{read_names[i]}', f'{j+1}', f'{reads[i][j]}', f'{qualities[i][j]}',f'{true_seq_vec[i][j]}' ])

    fh.close()

    return

nucs = ['A','C','G','T']
h = 1
reads = []
qualities = []
read_names = []
for record in SeqIO.parse("test.fastq", "fastq"):
    reads.append(str(record.seq))
    read_names.append(str(record.id))
    qualities.append(record.letter_annotations["phred_quality"])
print("Data read in.")
# Convert to Numpy Arrays
reads = np.asarray(reads)
qualities = np.asarray(qualities)
list_of_hstrings = []
sample_space = generate_hstring(list_of_hstrings, nucs, "", len(nucs), h)


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

pi = np.random.rand(4,4)
pi /= np.sum(pi, axis=0)
pi = pd.DataFrame(pi, index = ['A', 'C', 'G', 'T'], columns = ['A', 'C', 'G', 'T'])

# pi = np.array([[0.307758,  0.307826,  0.307808,  0.307868],
#                [0.181251,  0.181192,  0.181232,  0.181218],
#                [0.172631,  0.172585,  0.172591,  0.172574],
#                [0.338360,  0.338397,  0.338369,  0.338339]])
#
# pi = pd.DataFrame(pi, index = ['A', 'C', 'G', 'T'], columns = ['A', 'C', 'G', 'T'])

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
# # print("Tau", tau)
# print("Tau[z']", tau['A'])


# WHEN H= 1
# eta0 = np.array([0.00000000e+00, 7.20113984e-02, 0.00000000e+00, 7.55090035e-05,
#                  2.75119175e-04, 3.79419307e-04, 4.74624299e-04, 4.68815024e-04,
#                  7.44726913e-04, 6.95694663e-04, 4.81980457e-04, 5.74447950e-04,
#                  7.49614830e-04, 2.05038360e-03, 2.39501788e-03, 1.38230375e-03,
#                  1.31819380e-03, 1.74439781e-03, 1.91047621e-03, 2.33490108e-03,
#                  3.32658090e-03, 2.85067216e-03, 4.43019548e-03, 4.36144234e-03,
#                  6.18840510e-03, 6.21866901e-03, 8.33570387e-03, 1.12429593e-02,
#                  2.01422706e-02, 2.15923900e-02, 4.10163976e-02, 5.36052677e-02,
#                  1.27243980e-01, 5.38495418e-01, 6.08826244e-02, 0.00000000e+00,
#                  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00])
# eta1 = np.array([0.00000000e+00, 7.50862492e-02, 0.00000000e+00, 7.06168323e-05,
#                  2.98276566e-04, 3.87216443e-04, 4.96951248e-04, 4.58476698e-04,
#                  7.73873087e-04, 6.93415407e-04, 4.59388205e-04, 5.52726023e-04,
#                  7.43032256e-04, 2.04744842e-03, 2.44178890e-03, 1.41342701e-03,
#                  1.29793965e-03, 1.73162023e-03, 1.93691410e-03, 2.34223184e-03,
#                  3.33197319e-03, 2.96320536e-03, 4.48327944e-03, 4.40434677e-03,
#                  6.37230105e-03, 6.36685741e-03, 8.37792334e-03, 1.15871196e-02,
#                  2.03070198e-02, 2.20551400e-02, 4.11082603e-02, 5.38727261e-02,
#                  1.24279154e-01, 5.38382414e-01, 5.88766870e-02, 0.00000000e+00,
#                  0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00])

# pi = np.array([[0.307758,  0.307826,  0.307808,  0.307868],
#                [0.181251, 0.181192 , 0.181232,  0.181218],
#                [0.172631,  0.172585,  0.172591,  0.172574],
#                [0.338360,  0.338397,  0.338369,  0.338339]])

#
# pi = np.random.rand(4,4)
# pi /= np.sum(pi, axis=0)
# pi = pd.DataFrame(pi, index = ['A', 'C', 'G', 'T'], columns = ['A', 'C', 'G', 'T'])
#
# #
# tau = np.array([[0.350608,  0.377825,  0.046233,  0.225334],
#                 [0.248384,  0.033259,  0.220741,  0.497616],
#                 [0.212377,  0.166600,  0.257812,  0.363211],
#                 [0.245338,  0.280341,  0.099946,  0.374374]])
# tau = pd.DataFrame(tau, index = sample_space, columns = ['A', 'C', 'G', 'T'])
# gamma = np.array([0.26988402, 0.23379493, 0.13503349, 0.36128755])


# Verterbi(pi, gamma, tau, eta0, eta1, reads, qualities, h, sample_space, read_names)
# print("Tau", tau)


# for i in range(len(reads)):
#     for j in range(len(reads[0])):
#         output_string = ''
#         if reads[i][j] == true_seq[j]:
#             print("There are no errors")
#         else:
#             # print("Name of Read", i)
#             # print("Position of Error", j)
#             # print("Read base", reads[i][j])
#             # print("Quality Score of Base", qualities[i][j])
#             # print("Correct Base", true_seq[j])
#
#             output_string = "Name of Read: " + str(i) + " Position of Error: " + str(j) + " Read Base: " + reads[i][j] +
#             " Quality Score of Base: "+ str(qualities[i][j]) + " Correct Base: ", true_seq[j]
#










