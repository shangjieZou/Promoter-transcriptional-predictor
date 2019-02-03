from Try_new_vectoring_method_with_dict import main_convert
import numpy as np
from sklearn.cluster import KMeans
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from Build_Matrices import main_Matrices



dict_positive_numeric = main_convert()[0]
dict_constitutive_numeric = main_convert()[1]
dict_negative_numeric = main_convert()[2]
Total_seq_list = main_convert()[3]
edit_distance_matrix = main_Matrices()[0]
normalized_edit_dist_matrix = main_Matrices()[1]
Similarity_matrix_S = main_Matrices()[2]
Diagonal_matrix = main_Matrices()[3]
Laplacian_Matrix_G = main_Matrices()[4]


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
Eigen_Vector_reduced = get_eigen_val_vec(Laplacian_Matrix_G, cluster_nums)[1]


def export_eigvec():
    np.savetxt("Eigen_vector.csv", Eigen_Vector_reduced, delimiter=',')


# KMeans Clustering
k_means_cluster = KMeans(n_clusters=cluster_nums)
cluster_pred = k_means_cluster.fit_predict(Eigen_Vector_reduced)


# When hoping to get 2D fig, change 0 to 1
def two_dim_plot():
    Big_cluster = np.zeros((440, 2))
    Small_cluster = np.zeros((440, 2))
    j = 0
    for i in range(len(Eigen_Vector_reduced)):
        if cluster_pred[i] == 0:
            Big_cluster[j] = Eigen_Vector_reduced[i]
            j += 1
    for i in range(len(Eigen_Vector_reduced)):
        if cluster_pred[i] == 1:
            Small_cluster[j] = Eigen_Vector_reduced[i]
            j += 1
    plt.scatter(Big_cluster[:, 0], Big_cluster[:, 1], c="yellow", label = "cluster_0")
    plt.scatter(Small_cluster[:, 0], Small_cluster[:, 1], c="purple", label = "cluster_1")
    plt.xlim(-0.047673129462278, -0.0476731294622807)
    plt.legend()
    plt.show()


def three_dim_plot():
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(Eigen_Vector_reduced[:, 0], Eigen_Vector_reduced[:, 1],
               Eigen_Vector_reduced[:, 2], c=cluster_pred)
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
        if len(dict1[key]) > threshold_size * 4:
            print(key)
    for key in dict2.keys():
        if len(dict2[key]) > threshold_size * 4:
            print(key)
    for key in dict3.keys():
        if len(dict3[key]) > threshold_size * 4:
            print(key)



def statistics_on_clustering():
    list_index_for_category = []
    for i in range(len(Total_seq_list)):
        if i >= 125 and i <= 214:
            list_index_for_category.append("Constitutive")
        else:
            list_index_for_category.append("Inducible")
    Total_zero, Total_one, Total_constitutive_zero, Total_iducible_zero, \
    Total_constitutive_one, Total_iducible_one = 0, 0, 0, 0, 0, 0
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
    percentage_Consti_zero = (Total_constitutive_zero / Total_zero) * 100
    percentage_indusi_zero = (Total_iducible_zero / Total_zero) * 100
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


def get_seqs_in_cluster(cluster_value):
    cluster_x = []
    A = [1, 0, 0, 0]
    G = [0, 1, 0, 0]
    C = [0, 0, 1, 0]
    T = [0, 0, 0, 1]
    for i in range(len(cluster_pred)):
        if cluster_pred[i] == cluster_value:
            vectorized_seq = []
            for j in range(len(Total_seq_list[i])):
                if Total_seq_list[i][j] == 'a':
                    vectorized_seq += A
                elif Total_seq_list[i][j] == 'g':
                    vectorized_seq += G
                elif Total_seq_list[i][j] == 'c':
                    vectorized_seq += C
                elif Total_seq_list[i][j] == 't':
                    vectorized_seq += T
            cluster_x.append(vectorized_seq)
    return cluster_x


def get_ref_for_cluster(cluster_value):
    cluster_x = get_seqs_in_cluster(cluster_value)
    for i in range(len(cluster_x)):
        for key in dict_positive_numeric.keys():
            if dict_positive_numeric[key] == cluster_x[i]:
                print(key)
                print("positive")
        for key in dict_constitutive_numeric.keys():
            if dict_constitutive_numeric[key] == cluster_x[i]:
                print(key)
                print("constitutive")
        for key in dict_negative_numeric.keys():
            if dict_negative_numeric[key] == cluster_x[i]:
                print(key)
                print("negative")


def seq_len(seq_list):
    x = []
    for i in range(len(seq_list)):
        x.append(len(seq_list[i]) / 4)
    return x


def distribution_of_small_cluster():
    cluster_1 = get_seqs_in_cluster(1)
    len_distribute = seq_len(cluster_1)
    plt.hist(len_distribute, histtype='stepfilled', alpha=0.6, density=10, bins=40, color="red")
    plt.ylabel("Frequency/Class Interval")
    plt.xlabel("Length of Sequences")
    # plt.xlim((0,4000))
    plt.title("Distribution of Promoters' Length (cluster_1)")
    plt.show()


def distribution_of__big_cluster():
    cluster_0 = get_seqs_in_cluster(0)
    len_distribute = seq_len(cluster_0)
    plt.hist(len_distribute, histtype='stepfilled', alpha=0.6, density=10, bins=40, color="navy")
    plt.ylabel("Frequency/Class Interval")
    plt.xlabel("Length of Sequences")
    # plt.xlim((0,4000))
    plt.title("Distribution of Promoters' Length (cluster_0)")
    plt.show()



# try filtering large seqs
def cluster_after_filtering():
    def filtering_Eigvec():
        New_eigvec_list = []
        for i in range(len(Total_seq_list)):
            if len(Total_seq_list[i]) < 500:
                New_eigvec_list.append(Eigen_Vector_reduced[i])
        New_eigvec = np.array(New_eigvec_list)
        return New_eigvec

    New_eigvec_filtered = filtering_Eigvec()
    k_means_cluster_filtered = KMeans(n_clusters=cluster_nums)
    cluster_pred_filtered = k_means_cluster_filtered.fit_predict(New_eigvec_filtered)

    def two_dim_plot_filter():
        plt.scatter(New_eigvec_filtered[:, 0], New_eigvec_filtered[:, 1], c=cluster_pred_filtered)
        plt.xlim(-0.047673129462278, -0.0476731294622807)
        plt.show()

    def three_dim_plot_filtered():
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(New_eigvec_filtered[:, 0], New_eigvec_filtered[:, 1],
                   New_eigvec_filtered[:, 2], c=cluster_pred_filtered)
        plt.xlim(-0.047673129462278, -0.0476731294622807)
        ax.view_init(elev=10, azim=235)
        plt.show()

    two_dim_plot_filter()
    three_dim_plot_filtered()

