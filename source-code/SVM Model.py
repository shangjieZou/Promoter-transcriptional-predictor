from Try_new_vectoring_method_with_dict import main_convert
import numpy as np
from sklearn.model_selection import cross_val_score,train_test_split
from sklearn import svm
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score
from sklearn.model_selection import KFold
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


vectorized_positive_dict = main_convert()[0]
vectorized_constitutive_dict = main_convert()[1]
vectorized_negative_dict = main_convert()[2]

def filter_oversize(dict_seq, threshold):
    list_newdict_key = []
    list_newdict_value = []
    for key in dict_seq.keys():
        list_newdict_key.append(key)
        list_newdict_value.append(dict_seq[key])
    newdict = dict(zip(list_newdict_key, list_newdict_value))
    for key in dict_seq.keys():
        if len(dict_seq[key]) > threshold*4:
            del newdict[key]
    return newdict

reduced_vectorized_positive_dict = filter_oversize(vectorized_positive_dict, 500)  # len = 89
reduced_vectorized_constitutive_dict = filter_oversize(vectorized_constitutive_dict, 500) # len = 83
reduced_vectorized_negative_dict = filter_oversize(vectorized_negative_dict, 500)  # len = 203


# fill seq vectors with 0
def fill_zero(dict_seq):
    for key in dict_seq.keys():
        list_zero = [0]
        dict_seq[key] += (500*4-len(dict_seq[key]))*list_zero
    return dict_seq


# fill seqs to equal length
vectorized_positive_dict_final = fill_zero(reduced_vectorized_positive_dict)
vectorized_constitutive_dict_final = fill_zero(reduced_vectorized_constitutive_dict)
vectorized_negative_dict_final = fill_zero(reduced_vectorized_negative_dict)

vectorized_positive_list_final = [vectorized_positive_dict_final[key]
                                      for key in vectorized_positive_dict_final.keys()]
vectorized_constitutive_list_final = [vectorized_constitutive_dict_final[key]
                                      for key in vectorized_constitutive_dict_final.keys()]
vectorized_negative_list_final = [vectorized_negative_dict_final[key]
                                      for key in vectorized_negative_dict_final.keys()]


# transform to numpy array and divide training set as well as test set
Vector_whole = vectorized_constitutive_list_final + vectorized_positive_list_final + vectorized_negative_list_final
X = np.array(Vector_whole).reshape(375, 2000)
y_list_constitutive = 83*[0]
y_list_inducible = 292*[1]
# y_list_positive = 89*[1]
# y_list_negative = 203 * [2]

# y_list = y_list_positive+y_list_constitutive+y_list_negative
y_list = y_list_constitutive + y_list_inducible
y = np.array(y_list).reshape(375)


train_X,test_X,train_y,test_y = train_test_split(X,y,test_size=0.2,random_state=1)


def optimizing_C():
    C = 0.5
    list_C = []
    list_acc = []
    while C<48:
        predictive_model = svm.SVC(C=C, kernel='rbf', gamma='auto')
        predictive_model.fit(train_X, train_y)
        acc_cross = cross_val_score(predictive_model, test_X, test_y, cv = 7, scoring="accuracy")
        acc = sum(acc_cross)/7
        list_C.append(C)
        list_acc.append(acc)
        C += 0.5
    print(list_C)
    print(list_acc)
    max_acc = max(list_acc)
    dict_training = dict(zip(list_acc, list_C))
    max_C = dict_training[max_acc]
    print("max_C = "+str(max_C))
    print("max_acc = " + str(max_acc))

    plt.plot(list_C, list_acc, linewidth = '1', color='coral', linestyle=':', marker='|')
    plt.xlabel("Parameter C")
    plt.ylabel("Accuracy")
    plt.title("Optimizing Parameter C")
    plt.show()


def optimizing_gamma():
    C = 10
    gamma = 0.001
    list_gamma = []
    list_acc = []
    while gamma<1.0:
        predictive_model = svm.SVC(C=C, kernel='rbf', gamma= gamma)
        predictive_model.fit(train_X, train_y)
        prediction = predictive_model.predict(test_X)
        Accuracy = accuracy_score(prediction, test_y)
        list_gamma.append(gamma)
        list_acc.append(Accuracy)
        gamma += 0.001
    print(list_gamma)
    print(list_acc)
    max_acc = max(list_acc)
    dict_training = dict(zip(list_acc, list_gamma))
    max_gamma = dict_training[max_acc]
    print("max_gamma = "+str(max_gamma))
    print("max_acc = " + str(max_acc))

    plt.plot(list_gamma, list_acc, linewidth = '1', color='coral', linestyle=':', marker='|')
    plt.xlabel("Parameter gamma")
    plt.ylabel("Accuracy")
    plt.title("Optimizing Parameter gamma")
    plt.show()


def optimizing_gamma_C():
    C = 0
    list_gamma = []
    list_acc = []
    list_C = []
    while C<35:
        C += 0.5
        gamma = 0.001
        while gamma<0.5:
            predictive_model = svm.SVC(C=C, kernel='rbf', gamma= gamma)
            predictive_model.fit(train_X, train_y)
            prediction = predictive_model.predict(test_X)
            Accuracy = accuracy_score(prediction, test_y)
            list_gamma.append(gamma)
            list_C.append(C)
            list_acc.append(Accuracy)
            gamma += 0.001
            print("finished gamma = %s, C = %s"%(gamma, C))
    print(list_gamma)
    print(list_C)
    print(list_acc)
    max_acc = max(list_acc)
    dict_training_gamma = dict(zip(list_acc, list_gamma))
    dict_training_C = dict(zip(list_acc, list_C))
    max_gamma = dict_training_gamma[max_acc]
    max_C = dict_training_C[max_acc]
    print("max_gamma = "+str(max_gamma))
    print("max_C = "+str(max_C))
    print("max_acc = " + str(max_acc))


def try_svm_optimized_para():
    C = 10
    gamma = 0.015
    predictive_model = svm.SVC(C=C, kernel='rbf', gamma = gamma)
    acc_cross = cross_val_score(predictive_model, X, y, cv = 7, scoring="accuracy")
    print(acc_cross)
    acc = sum(acc_cross)/7
    print("acc = " + str(acc))


def Hold_out():
    C = 10
    gamma = 0.05
    predictive_model = svm.SVC(C=C, kernel='rbf', gamma=gamma)
    predictive_model.fit(train_X, train_y)
    prediction = predictive_model.predict(test_X)
    Accuracy = accuracy_score(prediction, test_y)
    print(Accuracy)  # 0.88


def cross_validation_K_fold():
    # k-fold
    C = 10
    gamma = 0.015
    for i in range(2,11):
        kf = KFold(n_splits=i, shuffle = True)
        predictive_model = svm.SVC(C=C, kernel='rbf', gamma=gamma)
        Accuracy_stack = 0
        for train, test in kf.split(X):
            predictive_model = predictive_model.fit(X[train], y[train])
            prediction = predictive_model.predict(X[test])
            Accuracy_by_module = accuracy_score(y[test],prediction)
            print(Accuracy_by_module)
            Accuracy_stack += Accuracy_by_module
        average_acc_by_module = Accuracy_stack/i
        print("%d fold cross validation"%(i))
        print("average_acc_by_module=" + str(average_acc_by_module))



"""
def plot_optimzing_three_dim():
    C = 0
    gamma = 0.05
    C = np.arange(0.5,48,0.5)
    list_gamma = []
    list_acc = []
    for i in range(len(C)):
        C_i = C[i]
        while gamma < 6:
            predictive_model = svm.SVC(C=C_i, kernel='rbf', gamma=gamma)
            predictive_model.fit(train_X, train_y)
            acc_cross = cross_val_score(predictive_model, test_X, test_y, cv=7, scoring='accuracy')
            acc = sum(acc_cross) / 7
            list_gamma.append(gamma)
            list_acc.append(acc)
            gamma += 0.05

    max_acc = max(list_acc)
    dict_training_gamma = dict(zip(list_acc, list_gamma))
    dict_training_C = dict(zip(list_acc, C))
    max_gamma = dict_training_gamma[max_acc]
    max_C = dict_training_C[max_acc]
    print("C:")
    print(C)
    print("list_gamma:")
    print(list_gamma)
    print("list_acc:")
    print(list_acc)
    print("max_acc = " + str(max_acc))
    print("max_gamma =a " + str(max_gamma))
    print("max_C = " + str(max_C))
    print(len(list_acc))
    print(len(list_gamma))
    print(len(C))
    Z = np.zeros((len(list_acc),3))
    for i in range(len(list_acc)):
        Z[i][0] = C[i]
        Z[i][1] = list_gamma[i]
        Z[i][2] = list_acc[i]
    Z.tofile("acc_across.dat", sep=',')
    print(Z)
"""
