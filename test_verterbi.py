
import numpy as np
from Bio import SeqIO
import pandas as pd

def Verterbi(pie, gamma, tau, eta0, eta1, reads, qualities):
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
        Zeta = np.zeros((len(reads[0]), 4))
        B = np.zeros((len(reads[0]),4))
        true_seq = ""
        for j in range(0, len(read)):
            # Forward Process
            for z in nucs:
                if j == 0:
                    # Intialize Zeta1(N) = gamma_N * eta0 or eta1 * pie
                    if read[j] == z:

                        Zeta[j][nucs.index(z)] = gamma[nucs.index(z)] * eta0[quality[j]-1] * pie[z][read[j]]

                    else:
                        Zeta[j][nucs.index(z)] = gamma[nucs.index(z)] * eta1[quality[j]-1] * pie[z][read[j]]
                else:
                    # Calculate new Zetas = eta0 or eta1 * pie * max(Zeta1(z') * tau(z', z)
                    old_z = code[np.argmax(Zeta[j-1])]
                    if read[j] == z:
                        Zeta[j][nucs.index(z)] = eta0[quality[j] - 1] * pie[z][read[j]] * np.max(Zeta[j-1] * tau[z][old_z])

                    else:
                        Zeta[j][nucs.index(z)] = eta1[quality[j] - 1] * pie[z][read[j]] * np.max(Zeta[j-1] * tau[z][old_z])
                    # Calculate a B matrix: Bj(z') = argmax(zeta_j-1(z') * tau(z', z) for j > 1
                    # print(Zeta[j-1] * tau[z][old_z])
                    B[j][nucs.index(z)] = np.argmax(Zeta[j-1] * tau[z][old_z])

        # for j in range(len(read)-1, 0, -1):
        #     # Backward Process:
        #     if j == len(read)-1:
        #         true_seq += np.argmax(Zeta[j])
        #     elif j > 1 and j < len(read)-1:
        #         # Zj = Bj+1(Zj+1)
        #         true_seq += B[j+1]
        print("Zeta:", Zeta)
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

tau = pd.DataFrame(np.array([[0.5, 0.25, 0, 0.25], [0, 0, 0.5, 0.5], [0.25, 0.25, 0.25, 0.25], [0.2, 0.3, 0.25, 0.25]]), index=['A', 'C', 'G', 'T'], columns=['A', 'C', 'G', 'T'])
gamma = np.array([0.1,0.2,0.4,0.3])
Verterbi(pi, gamma, tau, eta0, eta1, reads, qualities)