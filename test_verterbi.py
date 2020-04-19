
import numpy as np
from Bio import SeqIO
import pandas as pd

def Verterbi(pie, gamma, tau, eta0, eta1, reads, qualities, h, sample_space):




    for i in range(len(reads)):

    # Forward Part:
        read = reads[i]
        print("Read", read)
        quality = qualities[i]
        Zeta = np.ones((2, pow(4,h)))
        # B = np.zeros((len(reads[0]),pow(4,h)))
        B = np.empty((len(reads[0]), pow(4, h)), dtype= "str")
        true_seq = ""
        # for j in range(0, len(read)):
        # Forward Process
            # if j == 0:
        # for z in nucs:

        for z in sample_space:
            # Intialize Zeta1(N) = gamma_N * eta0 or eta1 * pie
            if read[0:h] == z:
                for l in range(0,h):

                    Zeta[0][sample_space.index(z)] *= gamma[nucs.index(z[l])] * eta0[quality[0]-1] * pie[z[l]][read[0]]

            else:
                Zeta[0][nucs.index(z)] = 0

        for j in range(h, len(read)):
                # Calculate new Zetas = eta0 or eta1 * pie * max(Zeta1(z') * tau(z', z)
                # print("J", j)
                # print("Read at position j", read[j])
                for z in sample_space:
                # for z in sample_space:

                    if read[j] == z[h-1]:
                        Zeta[1][sample_space.index(z)] = eta0[quality[j] - 1] * pie[z][read[j]] * np.max((Zeta[0] * tau[z]).to_numpy())

                    else:
                        Zeta[1][sample_space.index(z)] = eta1[quality[j] - 1] * pie[z][read[j]] * np.max((Zeta[0] * tau[z]).to_numpy())

                    # Calculate a B matrix: Bj(z') = argmax(zeta_j-1(z') * tau(z', z) for j > 1
                    B[j][sample_space.index(z)] = sample_space[np.argmax((Zeta[0] * tau[z]).to_numpy())]

                Zeta[0] = Zeta[1]

                # print(np.argmax(Zeta[j - 1]))
                # B[j][nucs.index(z)] = np.argmax(Zeta[j - 1])

        # print("B:", B)
        for j in range(len(read)-1, -1, -1):
            # Backward Process:
            if j == len(read)-1:
                # for z in sample_space:
                Z_hat = sample_space[np.argmax((Zeta[0] * tau[z]).to_numpy())]

                true_seq += Z_hat

            else:
                # Zj = Bj+1(Zj+1)
                new_Zhat = B[j+1][nucs.index(Z_hat)]

                true_seq += new_Zhat
                Z_hat = new_Zhat


        # Reverse the sequence
        print("Read seq", read)
        print("True seq", true_seq[::-1])

        # print("Zeta:", Zeta)
        # print("ARGmax:", np.argmax(Zeta, axis=1))
        # print("B:", B)
    return


reads = []
qualities = []
for record in SeqIO.parse("test.fastq", "fastq"):

    reads.append(str(record.seq))
    qualities.append(record.letter_annotations["phred_quality"])
print("Data read in.")
# Convert to Numpy Arrays
reads = np.asarray(reads)
qualities = np.asarray(qualities)


# eta0 = np.zeros(40)
# eta0[36:40] = [0.01, 0.19, 0.3, 0.5]
#
# eta1 = eta0 #np.random.dirichlet(np.ones(40))
# pi = np.array([[0.5, 0.25, 0.125, 0.125],[0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25]])
# pi = pd.DataFrame(pi, index= ('A', 'C', 'G', 'T'), columns= ('A', 'C', 'G', 'T'))
#
# tau = pd.DataFrame(np.array([[0.25, 0.25, 0.25, 0.25], [0.25, 0.1, 0.05, 0.05], [.03, 0.25, 0.47, 0.25], [0.1, 0.4, 0.25, 0.25]]), index=['A', 'C', 'G', 'T'], columns=['A', 'C', 'G', 'T'])
# gamma = np.array([0.1, 0.2, 0.4, 0.3])


eta0 = np.array([0.00000000e+00, 2.06474615e-03, 0.00000000e+00, 3.54430565e-06,
                 1.39900355e-04,  1.13642366e-04,  2.04139506e-04,  3.02724323e-04,
                 1.78070090e-04, 1.19373127e-04, 1.53961832e-04, 1.46672980e-04,
                 2.58022738e-04, 1.30748786e-03, 2.24574435e-03, 3.63709987e-04,
                 6.51856724e-04, 4.22038073e-04, 8.27234568e-04, 7.58945440e-04,
                 1.92969862e-03, 1.05492230e-03, 2.18821182e-03, 2.85799972e-03,
                 2.78445296e-03, 7.82761076e-04, 2.57061659e-03, 1.60179139e-03,
                 3.91137706e-03, 8.80759440e-03, 4.46991908e-02, 1.85757976e-02,
                 7.97212054e-02, 7.93891790e-01, 2.43607753e-02, 0.00000000e+00,
                 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00])
eta1 = np.array([0.00000000e+00, 1.36410686e-01, 0.00000000e+00, 1.30836536e-04,
                4.22968262e-04, 6.18844935e-04, 7.37817181e-04, 5.97800412e-04,
                1.27219026e-03, 1.18884281e-03, 7.33767859e-04, 9.13311192e-04,
                1.16396825e-03, 2.68607022e-03, 2.58700694e-03, 2.30158898e-03,
                1.86446193e-03, 2.86567696e-03, 2.87906942e-03, 3.70182084e-03,
                4.53667983e-03, 4.54958049e-03, 6.43267382e-03, 5.71419008e-03,
                9.36887556e-03, 1.11003044e-02, 1.33570034e-02, 2.00114139e-02,
                3.43413771e-02, 3.32285579e-02, 3.79699717e-02, 8.41302055e-02,
                1.64140798e-01, 3.18433534e-01, 8.96081045e-02, 0.00000000e+00,
                0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00])


pi = np.array([[0.580720 , 0.199521 , 0.303183 , 0.130883],
               [0.144903 , 0.298296 , 0.198463 , 0.060835],
               [0.125947 , 0.240025 , 0.261317 , 0.036037],
               [0.148430 , 0.262158 , 0.237037,  0.772245]])

pi = pd.DataFrame(pi, index = ['A', 'C', 'G', 'T'], columns = ['A', 'C', 'G', 'T'])

tau = np.array([[0.520965,  0.018157,  0.236503,  0.224375],
                [0.197305,  0.409745,  0.330057,  0.062894],
                [0.177568,  0.413912,  0.339237,  0.069282],
                [0.096983,  0.180927,  0.144058,  0.578033]])
tau = pd.DataFrame(tau, index = ['A', 'C', 'G', 'T'], columns = ['A', 'C', 'G', 'T'])
gamma = np.array([0.28861617,0.14515446,0.16969953, 0.39652983])

h = 1
nucs = ['A', 'C', 'G', 'T']
Verterbi(pi, gamma, tau, eta0, eta1, reads, qualities, h, nucs)
# print("Tau", tau)
# print("Tau[z']", tau['A'])














