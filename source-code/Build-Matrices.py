import numpy as np
from Try_new_vectoring_method_with_dict import main_convert
import math
from sklearn.cluster import KMeans
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt




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


#matrix_G(Diagonal_matrix, Similarity_matrix_S)
Laplacian_Matrix_G = np.fromfile("Laplacian_Matrix_G.dat", sep=',').reshape((len(Diagonal_matrix),
                                                                             len(Diagonal_matrix)))

def export_Laplacian_Matrix_G():
    print(Laplacian_Matrix_G)
    np.savetxt("Laplacian_Matrix_G.csv", Laplacian_Matrix_G, delimiter=',')


def get_eigen_val_vec(Matrix, cluster_nums):
    eigen_val, eigen_vec = np.linalg.eig(Matrix)
    dict_Eigenval = dict(zip(eigen_val, range(0, len(eigen_val))))
    sorted_eigval = np.sort(eigen_val)
    Lowest_k_eigval = []
    for i in range(cluster_nums):
        Lowest_k_eigval.append(sorted_eigval[i])
    index = []
    for key in Lowest_k_eigval:
        index.append(dict_Eigenval[key])
    Eigen_vec_selected = eigen_vec[:, index]
    return eigen_val[index], Eigen_vec_selected

# cluster_nums = 3
# cluster_nums = 4
cluster_nums = 2
Eigen_Vector_reduced = get_eigen_val_vec(Laplacian_Matrix_G,cluster_nums)[1]
def export_eigvec():
    np.savetxt("Eigen_vector.csv", Eigen_Vector_reduced, delimiter=',')


k_means_cluster = KMeans(n_clusters = cluster_nums)
cluster_pred = k_means_cluster.fit_predict(Eigen_Vector_reduced)


# When hoping to get 2D fig, change 0 to 1
def two_dim_plot():
    plt.scatter(Eigen_Vector_reduced[:,0], Eigen_Vector_reduced[:,1], c = cluster_pred)
    plt.xlim(-0.047673129462278, -0.0476731294622807)
    plt.legend()
    plt.show()


def three_dim_plot():
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(Eigen_Vector_reduced[:,0],Eigen_Vector_reduced[:,1],
               Eigen_Vector_reduced[:, 2], c = cluster_pred)
    plt.xlim(-0.047673129462278, -0.0476731294622807)
    ax.view_init(elev=10, azim=235)
    plt.show()


def find_out_strange_point():
    print(cluster_pred)
    print(type(cluster_pred))
    np.savetxt("cluster_pred.csv", cluster_pred, delimiter=',')





# Filter oversize seqs
def find_oversize(dict1, dict2, dict3, threshold_size):
    for key in dict1.keys():
        if len(dict1[key])> threshold_size * 4:
            print(key)
    for key in dict2.keys():
        if len(dict2[key]) > threshold_size * 4:
            print(key)
    for key in dict3.keys():
        if len(dict3[key]) > threshold_size * 4:
            print(key)


list_index_for_category = []
for i in range(len(Total_seq_list)):
    if i>=125 and i<=214:
        list_index_for_category.append("Constitutive")
    else:
        list_index_for_category.append("Inducible")


def statistics_on_clustering():
    Total_zero, Total_one, Total_constitutive_zero, Total_iducible_zero,\
    Total_constitutive_one, Total_iducible_one= 0,0,0,0,0,0
    for i in range(len(list_index_for_category)):
        if cluster_pred[i] == 0:
            Total_zero += 1
            if list_index_for_category[i] == "Constitutive":
                Total_constitutive_zero += 1
            else:
                Total_iducible_zero += 1

        else:
            Total_one += 1
            if list_index_for_category[i] == "Constitutive":
                Total_constitutive_one += 1
            else:
                Total_iducible_one += 1
    percentage_Consti_zero = (Total_constitutive_zero/Total_zero)*100
    percentage_indusi_zero = (Total_iducible_zero/Total_zero)*100
    percentage_Consti_one = (Total_constitutive_one / Total_one) * 100
    percentage_indusi_one = (Total_iducible_one / Total_one) * 100
    print("The Total number of Constitutive promoter marked as 0 is " + str(Total_constitutive_zero))
    print("The Total number of Indusible promoter marked as 0 is " + str(Total_iducible_zero))
    print("percentage of Constitutive promoter marked as 0 in all 0= %s percent" % (percentage_Consti_zero))
    print("percentage of Inducible promoter marked as 0 in all 0= %s percent" % (percentage_indusi_zero))
    print("The Total number of Constitutive promoter marked as 1 is " + str(Total_constitutive_one))
    print("The Total number of Indusible promoter marked as 1 is " + str(Total_iducible_one))
    print("percentage of Constitutive promoter marked as 1 in all 1= %s percent" % (percentage_Consti_one))
    print("percentage of Inducible promoter marked as 1 in all 1= %s percent" % (percentage_indusi_one))


def get_seqs_in_smaller_cluster():
    cluster_1 = []
    A = [1, 0, 0, 0]
    G = [0, 1, 0, 0]
    C = [0, 0, 1, 0]
    T = [0, 0, 0, 1]
    for i in range(len(cluster_pred)):
        if cluster_pred[i] == 1:
            vectorized_seq = []
            for j in range(len(Total_seq_list[i])):
                if Total_seq_list[i][j]== 'a':
                    vectorized_seq += A
                elif Total_seq_list[i][j] == 'g':
                    vectorized_seq += G
                elif Total_seq_list[i][j] == 'c':
                    vectorized_seq += C
                elif Total_seq_list[i][j] == 't':
                    vectorized_seq += T
            cluster_1.append(vectorized_seq)
    return cluster_1


def get_ref_for_cluster():
    cluster_1 = get_seqs_in_smaller_cluster()
    for i in range(len(cluster_1)):
        for key in dict_positive_numeric.keys():
            if dict_positive_numeric[key] == cluster_1[i]:
                print(key)
                print("positive")
        for key in dict_constitutive_numeric.keys():
            if dict_constitutive_numeric[key] == cluster_1[i]:
                print(key)
                print("constitutive")
        for key in dict_negative_numeric.keys():
            if dict_negative_numeric[key] == cluster_1[i]:
                print(key)
                print("negative")



def seq_len(seq_list):
    x = []
    for i in range(len(seq_list)):
        x.append(len(seq_list[i])/4)
    return x
cluster_1 = get_seqs_in_smaller_cluster()
len_distribute = seq_len(cluster_1)
plt.hist(len_distribute ,histtype='stepfilled', alpha=0.6, density=10, bins=40, color="red", label="Positive-regulated promoter")
plt.ylabel("Frequency/Class Interval")
plt.xlabel("Length of Sequences")
#plt.xlim((0,4000))
plt.title("Distribution of Promoters' Length (Smaller Cluster)")
plt.show()
