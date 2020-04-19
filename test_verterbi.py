
import numpy as np
from Bio import SeqIO
import pandas as pd

def Verterbi(pie, gamma, tau, eta0, eta1, reads, qualities, h):
    code = {0: "A", 1: "C", 2: "G", 3: "T"}
    # read = reads[0]
    # quality = qualities[0]
    nucs = ['A', 'C', 'G', 'T']
    # Zeta = np.zeros((2, 4))
    B = np.array(len(reads[0]))
    for i in range(len(reads)):
    # Forward Part:
        read = reads[i]
        quality = qualities[i]
        Zeta = np.zeros((2, pow(4,h)))
        B = np.zeros((len(reads[0]),pow(4,h)))
        true_seq = ""
        for j in range(0, len(read)):
            # Forward Process

            if j == 0:
                for z in nucs:
                    # Intialize Zeta1(N) = gamma_N * eta0 or eta1 * pie
                    if read[j] == z:
                        # print("Eta0:", gamma[nucs.index(z)] * eta0[quality[j]-1] * pie[z][read[j]])
                        # print("BEfore:", Zeta)
                        Zeta[0][nucs.index(z)] = gamma[nucs.index(z)] * eta0[quality[j]-1] * pie[z][read[j]]
                        # print("After",Zeta)

                    else:
                        # print("Eta1:", gamma[nucs.index(z)] * eta1[quality[j]-1] * pie[z][read[j]])
                        # print("Before:,", Zeta)
                        Zeta[0][nucs.index(z)] = gamma[nucs.index(z)] * eta1[quality[j]-1] * pie[z][read[j]]
                        # print("After", Zeta)
            else:
                # Calculate new Zetas = eta0 or eta1 * pie * max(Zeta1(z') * tau(z', z)
                for z in nucs:
                    old_z = code[np.argmax(Zeta[0])]
                    # print("Tau", tau)
                    # print("Tau(z')", tau[z])
                    if read[j] == z:

                        # print("After multiplying:", (Zeta[0] * tau[z]).to_numpy())
                        # print(np.max((Zeta[0] * tau[z]).to_numpy()))
                        Zeta[1][nucs.index(z)] = eta0[quality[j] - 1] * pie[z][read[j]] * np.max((Zeta[0] * tau[z]).to_numpy())

                    else:
                        # print("Max", np.max(Zeta[0] * tau[z][old_z]))
                        Zeta[1][nucs.index(z)] = eta1[quality[j] - 1] * pie[z][read[j]] * np.max((Zeta[0] * tau[z]).to_numpy())
                        # print("After multiplying:", (Zeta[0] * tau[z]).to_numpy())
                        # print(np.max((Zeta[0] * tau[z]).to_numpy()))

                    # Calculate a B matrix: Bj(z') = argmax(zeta_j-1(z') * tau(z', z) for j > 1
                    # B[j][nucs.index(z)] = np.argmax(Zeta[0] * tau[z][old_z])

                    B[j][nucs.index(z)] = np.argmax((Zeta[0] * tau[z]).to_numpy())


                Zeta[0] = Zeta[1]
                # print("B:", B)



                # print(np.argmax(Zeta[j - 1]))
                # B[j][nucs.index(z)] = np.argmax(Zeta[j - 1])


        # for j in range(len(read)-1, 0, -1):
        #     # Backward Process:
        #     if j == len(read)-1:
        #         true_seq += np.argmax(Zeta[j])
        #     elif j > 1 and j < len(read)-1:
        #         # Zj = Bj+1(Zj+1)
        #         true_seq += B[j+1]


        # print("Zeta:", Zeta)
        # print("ARGmax:", np.argmax(Zeta, axis=1))
        print("B:", B)
    return


reads = []
qualities = []
for record in SeqIO.parse("robs_test.fastq", "fastq"):

    reads.append(str(record.seq))
    qualities.append(record.letter_annotations["phred_quality"])
print("Data read in.")
# Convert to Numpy Arrays
reads = np.asarray(reads)
qualities = np.asarray(qualities)


eta0 = np.zeros(40)
eta0[36:40] = [0.01, 0.19, 0.3, 0.5]

eta1 = eta0 #np.random.dirichlet(np.ones(40))
pi = np.array([[0.5, 0.25, 0.125, 0.125],[0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25]])
pi = pd.DataFrame(pi, index= ('A', 'C', 'G', 'T'), columns= ('A', 'C', 'G', 'T'))

tau = pd.DataFrame(np.array([[0.25, 0.25, 0.25, 0.25], [0.25, 0.1, 0.05, 0.05], [.03, 0.25, 0.47, 0.25], [0.1, 0.4, 0.25, 0.25]]), index=['A', 'C', 'G', 'T'], columns=['A', 'C', 'G', 'T'])
gamma = np.array([0,0.3,0.4,0.3])
h = 1
Verterbi(pi, gamma, tau, eta0, eta1, reads, qualities, h)
# print("Tau", tau)
# print("Tau[z']", tau['A'])
