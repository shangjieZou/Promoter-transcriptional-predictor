import numpy as np
from Try_new_vectoring_method_with_dict import main_convert
import math


dict_positive_numeric = main_convert()[0]
dict_constitutive_numeric = main_convert()[1]
dict_negative_numeric = main_convert()[2]
Total_seq_list = main_convert()[3]


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
                                                   min(edit_distance_between_seqs[i - 1][j] + 1,
                                                       edit_distance_between_seqs[i][j - 1] + 1))
    return edit_distance_between_seqs[len_seq1][len_seq2]


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
"""
a = np.arange(8).reshape(4, 2)
a.tofile("edit_dinstance_matrix.dat", sep=',', format='%d')
b = np.fromfile("edit_dinstance_matrix.dat",dtype = np.int, sep = ',').reshape(4,2)
print(b)
"""

#distance_matrix()
edit_distance_matrix = np.fromfile("edit_dinstance_matrix.dat", dtype=np.int,
                                   sep=',').reshape((440,440))

def export_edit_dist_matrix():
    print(edit_distance_matrix)
    np.savetxt("edit_dinstance_matrix.csv", edit_distance_matrix, fmt='%d', delimiter=',')


def normalize_edit_dist_matrix(Matrix):
    M = np.zeros((len(Matrix), len(Matrix)))
    for i in range(len(Matrix)):
        for j in range(i, len(Matrix)):
            max_len_seq = max(len(Total_seq_list[i]), len(Total_seq_list[j]))
            print("max_len = "+str(max_len_seq))
            M[i][j] = Matrix[i][j]/max_len_seq
            print("M[i][j] = "+str(M[i][j]))
            M[j][i] = M[i][j]
            print("M[j][i] = " + str(M[j][i]))
    M.tofile("Normalized_edit_dist_matrix.dat", sep=',')

#normalize_edit_dist_matrix(edit_distance_matrix)

normalized_edit_dist_matrix = np.fromfile("Normalized_edit_dist_matrix.dat", sep=',').reshape((440,440))


def export_normalized_edit_dist():
    print(normalized_edit_dist_matrix)
    np.savetxt("Normalized_edit_dist_matrix.csv", normalized_edit_dist_matrix, delimiter=',')



# Similarity matrix built by Gaussian kernel
def matrix_S(Matrix):
    gamma = 0.3    # based on paper
    S = np.zeros((len(Matrix), len(Matrix)))
    M = normalized_edit_dist_matrix
    for i in range(len(Matrix)):
        for j in range(i, len(Matrix)):
            S[i][j] = math.exp(-1*M[i][j]/(2 * (gamma**2)))
            S[j][i] = S[i][j]
            print("S[i][j] = " + str(S[i][j]))
            print("S[j][i] = " + str(S[j][i]))
    S.tofile("Similarity_matrix_S.dat", sep=',')


#matrix_S(normalized_edit_dist_matrix)
Similarity_matrix_S = np.fromfile("Similarity_matrix_S.dat", sep=',').\
    reshape((len(normalized_edit_dist_matrix),len(normalized_edit_dist_matrix)))


def export_Similarity_matrix_S():
    print(Similarity_matrix_S)
    np.savetxt("Similarity_matrix_S.csv", Similarity_matrix_S, delimiter=',')



def Diagonal_matrix_D(Matrix):
    d = np.diag(np.zeros(len(Matrix)))
    for i in range(len(Matrix)):
        d[i][i] = sum(Matrix[i])
        print("d[%d][%d]=%s" % (i, i, d[i][i]))
    d.tofile("Diagonal_matrix_D.dat", sep=',')

#Diagonal_matrix_D(Similarity_matrix_S)

Diagonal_matrix = np.fromfile("Diagonal_matrix_D.dat", sep=',').reshape((440,440))

def export_Diagonal_matrix():
    print(Diagonal_matrix)
    np.savetxt("Diagonal_matrix.csv", Diagonal_matrix, delimiter=',')


# calculate eigendecomposition Matrix G (Laplacian Matrix)
def matrix_G(Diagonal_matrix, Similarity_matirx):
    G = np.zeros((len(Diagonal_matrix), len(Diagonal_matrix)))
    for i in range(len(Diagonal_matrix)):
        for j in range(len(Similarity_matirx)):
            G[i][j] = Diagonal_matrix[i][j] - Similarity_matirx[i][j]
            print("G[%d][%d]=%s" % (i, j, G[i][j]))
    G.tofile("Laplacian_Matrix_G.dat", sep=',')

