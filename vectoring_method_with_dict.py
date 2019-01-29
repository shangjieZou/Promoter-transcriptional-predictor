import matplotlib.pyplot as plt
import numpy as np


# This block is used to transform .txt into dicts
def convert_txt_to_dict(String_seq):
    BBa_name = []
    pure_seq = []
    String_seq_list = String_seq.split('\n')
    for i in range(len(String_seq_list)):
        list_one_line = String_seq_list[i]
        list_one_line_split = list_one_line.split(' : ')
        BBa_name.append(list_one_line_split[0])
        pure_seq.append(list_one_line_split[1])
    dict_seq = dict(zip(BBa_name, pure_seq))
    return dict_seq


def positive_seq_dict():
    f_positive_seq = open('Positive_seqs_with_ref_whole(125).txt', 'r')
    seq_positive = f_positive_seq.read()
    dict_seq_positive = convert_txt_to_dict(seq_positive)
    f_positive_seq.close()
    return dict_seq_positive

dict_seq_positive = positive_seq_dict()


def constitutive_seq_dict():
    f_constitutive_seq = open('Constitutive_seqs_with_ref_whole(90).txt', 'r')
    constitutive_seq = f_constitutive_seq.read()
    dict_seq_constitutive = convert_txt_to_dict(constitutive_seq)
    f_constitutive_seq.close()
    return dict_seq_constitutive

dict_seq_constitutive = constitutive_seq_dict()


def Negative_seq_dict():
    f_Negative_seq = open('Negative_seqs_with_ref_whole(226).txt', 'r')
    Negative_seq = f_Negative_seq.read()
    dict_seq_negative = convert_txt_to_dict(Negative_seq)
    f_Negative_seq.close()
    return dict_seq_negative

dict_seq_negative = Negative_seq_dict()




# This block is for screening seqs
def plot_distrbution_of_seq_len(dict_seq_positive, dict_seq_constitutive, dict_seq_negative):
    def accumulate_seq_len(dict_seq):
        x = []
        for value in dict_seq.values():
            x.append(len(value))
        return x
    x_positive = accumulate_seq_len(dict_seq_positive)
    x_constitutive = accumulate_seq_len(dict_seq_constitutive)
    x_negative = accumulate_seq_len(dict_seq_negative)
    plt.hist(x_positive,histtype='stepfilled', alpha=0.4, density=10, bins=40, color="navy", label="Positive-regulated promoter")
    plt.hist(x_constitutive, histtype='stepfilled', alpha=0.4, density=10, bins=40, color="green", label="Constitutive promoter")
    plt.hist(x_negative, histtype='stepfilled', alpha=0.4, density=10, bins=40, color="darkgoldenrod", label="Repressible promoter")
    plt.ylabel("Frequency/Class Interval")
    plt.xlabel("Length of Sequences")
    #plt.xlim((0,4000))
    plt.title("Distribution of Promoters' Length")
    plt.legend()
    plt.show()


def one_hot_vectorizing(dict_seq):
    A = [1,0,0,0]
    G = [0,1,0,0]
    C = [0,0,1,0]
    T = [0,0,0,1]
    list_key = []
    list_value = []
    for key in dict_seq.keys():
        list_key.append(key)
        vectorized_seq = []
        for i in range(len(dict_seq[key])):

            if dict_seq[key][i] == 'a':
                vectorized_seq += A
            elif dict_seq[key][i] == 'g':
                vectorized_seq += G
            elif dict_seq[key][i] == 'c':
                vectorized_seq += C
            elif dict_seq[key][i] == 't':
                vectorized_seq += T
        list_value.append(vectorized_seq)
    dict_vectorized = dict(zip(list_key, list_value))
    return dict_vectorized

dict_positive_numeric = one_hot_vectorizing(dict_seq_positive)
dict_constitutive_numeric = one_hot_vectorizing(dict_seq_constitutive)
dict_negative_numeric = one_hot_vectorizing(dict_seq_negative)

"""
for key in dict_positive_numeric.keys():
    print(key)
    print(dict_positive_numeric[key])
"""
"""
for key in dict_seq_negative.keys():
    print(key)
    print(dict_seq_negative[key])
"""



# This block is used to generate Distant Matrix
def edit_distance(seq1, seq2):
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)
    edit_distance_between_seqs = np.zeros((len_seq1 + 1, len_seq2 + 1))
    for i in range(len_seq1 + 1):
        edit_distance_between_seqs[i][0] = i
    for j in range(len_seq2 + 1):
        edit_distance_between_seqs[0][j] = j

    for i in range(1, len_seq1 + 1):
        for j in range(1, len_seq2 + 1):
            delta = 0 if seq1[i - 1] == seq2[j - 1] else 1
            edit_distance_between_seqs[i][j] = min(edit_distance_between_seqs[i - 1][j - 1] + delta,
                                                   min(edit_distance_between_seqs[i - 1][j] + 1, edit_distance_between_seqs[i][j - 1] + 1))
    return edit_distance_between_seqs[len_seq1][len_seq2]


Total_seq_list = []
for key in dict_seq_positive.keys():
    Total_seq_list.append(dict_seq_positive[key])
for key in dict_seq_constitutive.keys():
    Total_seq_list.append(dict_seq_constitutive[key])
for key in dict_seq_negative.keys():
    Total_seq_list.append(dict_seq_negative[key])
"""
print(len(Total_seq_list))
for i in range(len(Total_seq_list)):
    print(Total_seq_list[i])
    print('\n')
"""


def distance_matrix():
    dim = len(dict_positive_numeric) + len(dict_constitutive_numeric) + len(dict_negative_numeric)
    D = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(i, dim):
            D[i][j] = edit_distance(Total_seq_list[i], Total_seq_list[j])
            D[j][i] = D[i][j]
            if(1):
                print("D[%s][%d] finished"%(i, j))
                print(D[i][j])
                print("D[%s][%d]=%d"%(j, i,D[j][i]))
    D.tofile("edit_dinstance_matrix.dat", sep=',', format='%d')

distance_matrix()
"""
a = np.arange(8).reshape(4, 2)
a.tofile("edit_dinstance_matrix.dat", sep=',', format='%d')
b = np.fromfile("edit_dinstance_matrix.dat",dtype = np.int, sep = ',').reshape(4,2)
print(b)
"""