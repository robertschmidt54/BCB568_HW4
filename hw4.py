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
                # print(tau[s[h-1]][t])

                alpha_vec[j-h][sample_space.index(s)] += alpha_vec[j-h-1][sample_space.index(t)]*tau[s[h-1]][t]
            # read[j] == s[h-1]:
            if read[j-h] == s[h-1]:
                alpha_vec[j-h][sample_space.index(s)] *= eta0[quality[j-h]-1]* pi[s[h-1]][read[j-h]]
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
        # print("z", z)
        # print("Temp_e", temp_e)

        e_num[sample_space.index(z)] = temp_e

    sum_eijzw = np.sum(e_num)
    # print("E_num", e_num)
    if sum_eijzw != 0:
        e = e_num/sum_eijzw
    else:
        e = e_num
    # print("After Division:", e)
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

def simplifyE(E_ijn, samplespace):
    simple_E = np.zeros((4, len(reads), len(reads[0])))
    for s in samplespace:
        if s.endswith('A'):
            simple_E[0] += E_ijn[sample_space.index(s)]
        elif s.endswith('C'):
            simple_E[1] += E_ijn[sample_space.index(s)]
        elif s.endswith('G'):
            simple_E[2] += E_ijn[sample_space.index(s)]
        else:
            simple_E[3] += E_ijn[sample_space.index(s)]

    return simple_E

def Update_Gamma(ei1n, y, h):

    # temp_gamma = np.zeros(pow(4,h))
    # temp_gamma += ei1n
    return np.sum(ei1n, axis=1)/y



def Update_Pi(E_ijn, sample_space, reads, h):
    if h > 1:
        E_ijn = simplifyE(E_ijn, sample_space)

    # print("E_ijn", E_ijn)
    piAn = np.zeros(4)
    piCn = np.zeros(4)
    piGn = np.zeros(4)
    piTn = np.zeros(4)
    y = 0 #Counter for reads because numpy doesn't have a good index function.

    for read in reads:
        k = 0  # Counter for position in read because reasons.
        for j in read:
            if j == 'A':
                piAn += E_ijn[:,y,k] #[E_ijn[0, y, k], E_ijn[1, y, k], E_ijn[2,y,k], E_ijn[3,y,k]]
            elif j == 'C':
                piCn += E_ijn[:,y,k] #[E_ijn[0, y, k], E_ijn[1, y, k], E_ijn[2,y,k], E_ijn[3,y,k]]
            elif j == 'G':
                piGn += E_ijn[:,y,k] #[E_ijn[0, y, k], E_ijn[1, y, k], E_ijn[2,y,k], E_ijn[3,y,k]]
            else:
                piTn += E_ijn[:,y,k] #[E_ijn[0, y, k], E_ijn[1, y, k], E_ijn[2, y, k], E_ijn[3, y, k]]
            k += 1
        y += 1


    e_sums = np.sum(E_ijn, axis=(1,2))
    # print("E_sums: ", e_sums)
    # print("E_ijT: ", E_ijn[3,:,:])
    # print("E_sum of T: ", np.sum(E_ijn[3,:,:]))
    pi_temp = np.stack((piAn, piCn, piGn, piTn),axis=0)
    pi_temp /= e_sums

    pi_out = pd.DataFrame(pi_temp , index=['A', 'C', 'G', 'T'], columns = ['A', 'C', 'G', 'T'])
    return pi_out






def Update_Eta(reads, qualities, E_ijn, sample_space):
    '''
        loop through reads and posistions:
            if r == 'A':
                eta0_temp[quality -1] += E_ijn[A, read number, posistion number]
                eta1_temp[quality -1] += E_ijn[C, T, and G, read number, posistion number]
            Similar if statements for C, G, T.

        eta0_temp/sum(eta0_temp)
        eta1_temp/sum(eta1_temp)
        '''
    if h > 1:
        E_ijn =simplifyE(E_ijn, sample_space)
    eta0_temp = np.zeros(40)
    eta1_temp = np.zeros(40)
    for i in range(len(reads)):
        quality = qualities[i]
        read = reads[i]
        for j in range(len(read)):
            q = quality[j]
            r = read[j]

            if r == 'A':
                eta0_temp[q-1] += E_ijn[0, i, j]
                eta1_temp[q-1] += E_ijn[1, i, j] + E_ijn[2, i, j] + E_ijn[3, i, j]
            elif r == 'C':
                eta0_temp[q-1] += E_ijn[1, i, j]
                eta1_temp[q-1] += E_ijn[0, i, j] + E_ijn[2, i, j] + E_ijn[3, i, j]
            elif r == 'G':
                eta0_temp[q-1] += E_ijn[2, i, j]
                eta1_temp[q-1] += E_ijn[1, i, j] + E_ijn[0, i, j] + E_ijn[3, i, j]
            else:
                eta0_temp[q-1] += E_ijn[3, i, j]
                eta1_temp[q-1] += E_ijn[1, i, j] + E_ijn[2, i, j] + E_ijn[0, i, j]

    return eta0_temp/sum(eta0_temp), eta1_temp/sum(eta1_temp)


def Update_Tau(eijzw, sample_space):

    eijzw_sums = eijzw.sum(axis= 1, keepdims=True)
    new_tau = eijzw/eijzw_sums

    return pd.DataFrame(new_tau, index= sample_space, columns= ('A', 'C', 'G', 'T'))


def Convergence(new_gamma, old_gamma,new_pie, old_pie, new_eta0,old_eta0, new_eta1, old_eta1, new_tau, old_tau):

    global pie_is_con
    global eta_is_con
    global tau_is_con
    global gamma_is_con

    if max(new_gamma - old_gamma) <= (10**(-6)):
        print("")
        print("Max Gamma Difference", np.max(new_gamma - old_gamma))
        print("")
        gamma_is_con = True


    if (new_pie - old_pie).values.max() <= (10**(-3)):
        print("")
        print("Max Pie Difference", (new_pie - old_pie).values.max())
        print("")
        pie_is_con = True

    if max(new_eta1 - old_eta1) <= 10**(-3) and max(new_eta0 - old_eta0) <= (10**(-3)):
        print("")
        print("Eta1 Difference", max(new_eta1 - old_eta1))
        print("")
        print("")
        print("Eta0 Difference", max(new_eta0 - old_eta0))
        print("")
        eta_is_con = True

    if (new_tau - old_tau).values.max() <= 10**(-3):
        print("")
        print("Tau Difference", (new_tau - old_tau).values.max())
        print("")
        tau_is_con = True

    return gamma_is_con and pie_is_con and eta_is_con and tau_is_con



code={"A":0, "C":1, "G":2, "T":3 }
nucs = ['A', 'C', 'G', 'T']

# Import Fastq
reads = []
qualities = []
print("Reading in Data.")
for record in SeqIO.parse("errored_reads.fastq", "fastq"):

    reads.append(str(record.seq))
    qualities.append(record.letter_annotations["phred_quality"])
print("Data read in.")
# Convert to Numpy Arrays
reads = np.asarray(reads)
qualities = np.asarray(qualities)


# Da Real stuff
# Order of markov chain
h = 7
# Run EM algorithm
pie_is_con = False
eta_is_con = False
tau_is_con = False
gamma_is_con = False
# Initial state matrix
new_tau = intialize_Tau(nucs, h)

new_gamma = np.random.dirichlet(np.ones(pow(4,h)))
new_eta0 = np.random.dirichlet(np.ones(40))
new_eta1 = np.random.dirichlet(np.ones(40))
pi = np.random.rand(4,4)
pi /= np.sum(pi, axis=0)
new_pi = pd.DataFrame(pi, index = ['A', 'C', 'G', 'T'], columns = ['A', 'C', 'G', 'T'])
list_of_hstrings = []
sample_space = generate_hstring(list_of_hstrings, nucs, "", len(nucs), h)

old_tau = np.zeros((pow(4,h), 4))
old_tau = pd.DataFrame(old_tau, index = sample_space, columns = ['A', 'C', 'G', 'T'])
old_gamma = np.zeros(pow(4,h))
old_eta0 = np.zeros(40)
old_eta1 = np.zeros(40)
old_pi = np.zeros((4,4))
old_pi = pd.DataFrame(old_pi,  index = ['A', 'C', 'G', 'T'], columns = ['A', 'C', 'G', 'T'])


# Main Function
iteration = 0
while not Convergence(new_gamma, old_gamma, new_pi, old_pi, new_eta0, old_eta0, new_eta1, old_eta1, new_tau, old_tau):
    E_ijzw = np.zeros((pow(4, h), 4))
    E_ijn = np.ones((pow(4, h), len(reads), len(reads[0])))
    # print("Test COnvergence")
    # print("New_gamma:", max(new_gamma - old_gamma))
    # print("New_Pi", (new_pi - old_pi).values.max())
    # print("New_Tau", (new_tau - old_tau).values.max())
    # print("New_Eta0", max(new_eta0 - old_eta0))
    # print("New_Eta1", max(new_eta1 - old_eta1))
    print(iteration)

    if iteration > 10:
        break
    for i in range(0, len(reads)):
        # Calculate E Function
        if i % 100 == 0:
            print(str(i) + " Reads processed.")

        temp_eijn = e_func(new_eta0, new_eta1, new_pi, new_gamma, reads[i], qualities[i], new_tau, h)[0]
        E_ijzw += e_func(new_eta0, new_eta1, new_pi, new_gamma, reads[i], qualities[i], new_tau, h)[1]

        for z in sample_space:
            E_ijn[sample_space.index(z)][i] = temp_eijn[z]

        # Calculating the nums as we go

    # # M step Updates
    if not gamma_is_con:
        old_gamma = new_gamma
        new_gamma = Update_Gamma(E_ijn[:,:,0], len(reads), h)

    if not tau_is_con:
        old_tau = new_tau
        new_tau = Update_Tau(E_ijzw, sample_space)

    if not pie_is_con:
        old_pi = new_pi
        new_pi = Update_Pi(E_ijn, sample_space, reads, h)
    if not eta_is_con:
        old_eta1 = new_eta1
        old_eta0 = new_eta0
        new_eta0, new_eta1 = Update_Eta(reads, qualities, E_ijn, sample_space)
    iteration += 1


print("Final Gamma Estimate:\n", new_gamma)
print("Final Pi Estimate\n", new_pi)
print("Final Tau Estimate\n", new_tau)
print("Final Eta0 Estimate\n", new_eta0)
print("Final Eta1 Estimate\n", new_eta1)
print(iteration)

# Write to file for Verterbi
# h = 7
fh = open('Parameter_estimates', 'w')
fh.write("H-order: ")
fh.write(h)
fh.write("\n")
fh.write("Final Gamma Estimate:")
fh.write(new_gamma)
fh.write("\n")
fh.write("Final Pi Esimate:")
fh.write(new_pi)
fh.write("\n")
fh.write("Final Tau Estimate:")
fh.write(new_tau)
fh.write("\n")
fh.write("Final Eta0 Estimate")
fh.write(new_eta0)
fh.write("\n")
fh.write("Final Eta1 Estimate")
fh.write(new_eta1)

fh.close()



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



# Verterbi algorithm
# def Verterbi():

# Forward Part:
# Intialize Zeta1(N) = gamma_N * eta0 or eta1 * pie
# Calculate new Zetas = eta0 or eta1 * pie * max(Zeta1(z') * tau(z', z)
# Calculate a B matrix: Bj(z') = argmax(zeta_j-1(z') * tau(z', z) for j > 1

# Backward Part:
# Z_2 = argmax( B_2(A), B_2(C), B_2(G), B_2(T))
# Z_1 = B_2(Z_2)

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
        B = np.zeros((len(reads[0],4)))
        true_seq = ""
        for j in range(0, len(read)):
            # Forward Process
            for z in nucs:
                if j == 0:
                    # Intialize Zeta1(N) = gamma_N * eta0 or eta1 * pie
                    if read[j] == z:
                        Zeta[j][nucs.index(z)] = gamma[nucs.index(z)] * eta0[quality-1] * pie[z][read[j]]

                    else:
                        Zeta[j][nucs.index(z)] = gamma[nucs.index(z)] * eta1[quality-1] * pie[z][read[j]]
                else:
                    # Calculate new Zetas = eta0 or eta1 * pie * max(Zeta1(z') * tau(z', z)
                    old_z = code[np.argmax(Zeta[j-1])]
                    if read[j] == z:
                        Zeta[j][nucs.index(z)] = eta0[quality - 1] * pie[z][read[j]] * np.max(Zeta[j-1] * tau[z][old_z])

                    else:
                        Zeta[j][nucs.index(z)] = eta1[quality - 1] * pie[z][read[j]] * np.max(Zeta[j-1] * tau[z][old_z])
                    # Calculate a B matrix: Bj(z') = argmax(zeta_j-1(z') * tau(z', z) for j > 1
                    B[j][nucs.index(z)] = np.argmax(Zeta[j-1] * tau[z][old_z])

        for j in range(len(read)-1, 0, -1):
            # Backward Process:
            if j == len(read)-1:
                true_seq += np.argmax(Zeta[j])
            elif j > 1 and j < len(read)-1:
                # Zj = Bj+1(Zj+1)
                true_seq += B[j+1]
    pass





# Z_2 = argmax( B_2(A), B_2(C), B_2(G), B_2(T))
# Z_1 = B_2(Z_2)
