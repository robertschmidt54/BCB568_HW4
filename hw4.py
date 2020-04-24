#### Imports ####
import numpy as np
from Bio import SeqIO
import pandas as pd
import csv


#### Functions ####
def intialize_Tau(nucs, h):
    """
    This function makes a first guess at the tau estimate
    """
    list_of_hstrings = []
    generate_hstring(list_of_hstrings, nucs, "", len(nucs), h)

    matrix = np.random.rand(pow(4,h), 4)
    test = matrix / matrix.sum(axis=1)[:, None]

    Tau = test
    return Tau


def generate_hstring(h_list, set_char, prefix, set_len, str_len):
    """
    This function generates the sample space of the Hidden sequence.

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


def alpha_func(eta0, eta1, pi, gamma, read, quality, tau, h, sample_space):
    """
    Forward algorithm in the Baum Welch algorithm.
    """

    alpha_vec = np.zeros((len(read) - h, pow(4, h)))

    s = read[0:h]
    alpha_vec[0][sample_space.index(s)] = 1
    # Assume the first H-mer is error free
    for l in range(0, h):
        # Loop through first h-mer
        alpha_vec[0][sample_space.index(s)] *= eta0[quality[l]-1] * pi[s[l]][read[l]] * gamma[code[s[l]]]

    for j in range(h+1, len(read)):

        for s in sample_space:
            for t in sample_space:
                alpha_vec[j-h][sample_space.index(s)] += alpha_vec[j-h-1][sample_space.index(t)]*tau[s][t]

            if read[j-1] == s[h-1]:
                alpha_vec[j-h][sample_space.index(s)] *= eta0[quality[j-1]-1] * pi[s[h-1]][read[j-1]]

            else:
                alpha_vec[j-h][sample_space.index(s)] *= eta1[quality[j-1]-1] * pi[s[h-1]][read[j-1]]

    return alpha_vec


def eijzw_num(alpha, beta, pi, eta0, eta1, tau, r, q, sample_space):
    """
    Calculates the numerator of eijzw matrix.
    """
    e_num = np.zeros((pow(4,h),pow(4,h)))
    # zloop
    for z in sample_space:
        # wloop
        temp_e = np.ones(pow(4,h))

        for w in sample_space:
            # Equation to calculate num of eijzw: alpha(z) * beta(w) * eta0(r=w) or eta1(r != w) * tau zw * pi(r,w)
            if r == w[h-1]:
                temp_e[sample_space.index(w)] = alpha[sample_space.index(z)] * beta[sample_space.index(w)] * eta0[q-1] * tau[w][z] * pi[w[-1]][r]
            else:
                temp_e[sample_space.index(w)] = alpha[sample_space.index(z)] * beta[sample_space.index(w)] * eta1[q- 1] * tau[w][z] * pi[w[-1]][r]

        e_num[sample_space.index(z)] = temp_e

    sum_eijzw = np.sum(e_num)

    if sum_eijzw != 0:
        e = e_num/sum_eijzw
    else:
        e = e_num
    return e




def e_func(eta0, eta1, pi, gamma, read, quality, tau, h, sample_space):

    """
    Calculates Beta matrix.
    This function determines our eijn and eijzw for each read using the Baum Welch's
    forward (alpha) and backward(beta) algorithms.
    """

    beta_out = np.ones((2, pow(4, h)))
    alpha = alpha_func(eta0, eta1, pi, gamma, read, quality, tau, h, sample_space)

    E_ijn = dict.fromkeys(sample_space, [])

    E_ijzw = np.zeros((pow(4, h), pow(4,h)))
    e = np.zeros((len(read), pow(4, h)))
    for j in range(len(read) - 1, -1, -1):
        # Calculate beta

        for s in sample_space:
            beta_out[0][sample_space.index(s)] = 0
            for t in sample_space:

                if read[j] == s[h - 1]:
                    beta_out[0][sample_space.index(s)] += beta_out[1][sample_space.index(t)] * eta0[quality[j] - 1] * pi[t[h - 1]][read[j]] * tau[t][s]
                else:
                    beta_out[0][sample_space.index(s)] += beta_out[1][sample_space.index(t)] * eta1[quality[j] - 1] * pi[t[h - 1]][read[j]] * tau[t][s]

        # Calculate the numerator of E_ijn for 1 read
        e[j] = np.multiply(alpha[j-h], beta_out[1])


        # Calculate the numerator of E_ijzw for 1 read
        a = alpha[j-h-1]
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

        for b, m in zip(E_ijn[z], eijn_sums):
            if b != 0:

                E_ijn[z] = b/m

    return E_ijn, E_ijzw


def simplifyE(E_ijn, samplespace):
    """
    Condenses a 4^h by 4^h matrix down to a 4 by 4.
    """
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
    """
    Estimates a new Gamma vector
    """
    return np.sum(ei1n, axis=1)/y



def Update_Pi(E_ijn, sample_space, reads, h):
    """
    Estimates a new Pi matrix
    """
    if h > 1:
        E_ijn = simplifyE(E_ijn, sample_space)

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

    pi_temp = np.stack((piAn, piCn, piGn, piTn),axis=0)
    pi_temp = pi_temp / pi_temp.sum(axis=0)[None,:]


    pi_out = pd.DataFrame(pi_temp , index=['A', 'C', 'G', 'T'], columns = ['A', 'C', 'G', 'T'])
    return pi_out



def Update_Eta(reads, qualities, E_ijn, sample_space):
    """
    Estimates a new eta0 and eta1 vector.
    Pseudocode:
    loop through reads and posistions:
            if r == 'A':
                eta0_temp[quality -1] += E_ijn[A, read number, posistion number]
                eta1_temp[quality -1] += E_ijn[C, T, and G, read number, posistion number]
            Similar if statements for C, G, T.

        eta0_temp/sum(eta0_temp)
        eta1_temp/sum(eta1_temp)

    """

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
    """
    Estimates a new tau matrix.
    """

    eijzw_sums = eijzw.sum(axis= 1, keepdims=True)
    new_tau = eijzw/eijzw_sums

    return pd.DataFrame(new_tau, index= sample_space, columns= sample_space)


def Convergence(new_gamma, old_gamma,new_pie, old_pie, new_eta0,old_eta0, new_eta1, old_eta1, new_tau, old_tau):
    """
    Test the Convergence criteria after each iteration of EM algorithm.
    """
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


def Verterbi(pie, gamma, tau, eta0, eta1, reads, qualities, h, sample_space, read_names):
    """
    Determines the final true sequence of each read. Outputs the errors in a .csv file
    """

    true_seq_vec = []
    with open('corrected_error_final_data.csv', mode = 'w') as fh:
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




#### Main Function ####
## Intialize Variables ##
code={"A":0, "C":1, "G":2, "T":3 }
nucs = ['A', 'C', 'G', 'T']

# Import Fastq
reads = []
qualities = []
read_names = []
print("Reading in Data.")
for record in SeqIO.parse("test.fastq", "fastq"):
    read_names.append(str(record.id))
    reads.append(str(record.seq))
    qualities.append(record.letter_annotations["phred_quality"])
print("Data read in.")
# Convert to Numpy Arrays
reads = np.asarray(reads)
qualities = np.asarray(qualities)

# Order of markov chain
h = 2

# Convergence variables
pie_is_con = False
eta_is_con = False
tau_is_con = False
gamma_is_con = False

# generate sample space
list_of_hstrings = []
sample_space = generate_hstring(list_of_hstrings, nucs, "", len(nucs), h)

# Parameters
# Initial state matrix
new_gamma = np.random.dirichlet(np.ones(pow(4,h)))

# Emission distribution
new_eta0 = np.random.dirichlet(np.ones(40))
new_eta1 = np.random.dirichlet(np.ones(40))
pi = np.random.rand(4,4)
pi = pi/pi.sum(axis=0)[None,:]
new_pi = pd.DataFrame(pi, index = ['A', 'C', 'G', 'T'], columns = ['A', 'C', 'G', 'T'])


# Transition probability
tau = intialize_Tau(nucs, h)
new_tau = np.zeros((pow(4, h), pow(4, h)))
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
new_tau = pd.DataFrame(new_tau, index = sample_space, columns = sample_space)

old_tau = np.zeros((pow(4,h), pow(4,h)))
old_tau = pd.DataFrame(old_tau, index = sample_space, columns = sample_space)
old_gamma = np.zeros(pow(4,h))
old_eta0 = np.zeros(40)
old_eta1 = np.zeros(40)
old_pi = np.zeros((4,4))
old_pi = pd.DataFrame(old_pi,  index = ['A', 'C', 'G', 'T'], columns = ['A', 'C', 'G', 'T'])


iteration = 0
## Run EM algorithm ##
while not Convergence(new_gamma, old_gamma, new_pi, old_pi, new_eta0, old_eta0, new_eta1, old_eta1, new_tau, old_tau):
    E_ijzw = np.zeros((pow(4, h), pow(4,h)))
    E_ijn = np.ones((pow(4, h), len(reads), len(reads[0])))

    print(iteration)

    if iteration > 10:
        break
    for i in range(0, len(reads)):
        # Calculate E Function
        if i % 100 == 0:
            print(str(i) + " Reads processed.")

        temp_eijn = e_func(new_eta0, new_eta1, new_pi, new_gamma, reads[i], qualities[i], new_tau, h, sample_space)[0]
        E_ijzw += e_func(new_eta0, new_eta1, new_pi, new_gamma, reads[i], qualities[i], new_tau, h, sample_space)[1]

        for z in sample_space:
            E_ijn[sample_space.index(z)][i] = temp_eijn[z]

    # M step Updates
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

## Output Estimates and Errors found in the data ##
print("Final Gamma Estimate:\n", new_gamma)
print("Final Pi Estimate\n", new_pi)
print("Final Tau Estimate\n", new_tau)
print("Final Eta0 Estimate\n", new_eta0)
print("Final Eta1 Estimate\n", new_eta1)
print(iteration)

# Run Verterbi Algorithm to Output True Sequences
Verterbi(new_pi, new_gamma, new_tau, new_eta0, new_eta1, reads, qualities, h, sample_space, read_names)

