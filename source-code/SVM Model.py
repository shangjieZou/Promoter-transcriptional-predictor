from VectorizingSequences import Vectoring_main
import numpy as np
from sklearn.model_selection import cross_val_score,train_test_split
from sklearn import svm
import matplotlib.pyplot as plt

vectorized_positive_list = Vectoring_main()[0]
vectorized_constitutive_list = Vectoring_main()[1]
vectorized_negative_list = Vectoring_main()[2]


# fill seq vectors with 0
def fill_zero(list_seq):
    for i in range(len(list_seq)):
        list_zero = [0]
        list_seq[i] += (800-len(list_seq[i]))*list_zero
    return list_seq

vectorized_positive_list_final = fill_zero(vectorized_positive_list)
vectorized_constitutive_list_final = fill_zero(vectorized_constitutive_list)
vectorized_negative_list_final = fill_zero(vectorized_negative_list)


# transform to numpy array and divide training set as well as test set
Vector_whole = vectorized_constitutive_list_final + vectorized_positive_list_final + vectorized_negative_list_final
X = np.array(Vector_whole).reshape(201, 800)
y_list_constitutive = 42*[0]
y_list_regulated = 159*[1]
#y_list_positive = 57*[1]
#y_list_negative = 102 * [2]

#y_list = y_list_positive+y_list_constitutive+y_list_negative
y_list = y_list_constitutive + y_list_regulated
y = np.array(y_list).reshape(201)

def normal_training():
    train_X,test_X,train_y,test_y = train_test_split(X,y,test_size=0.2,random_state=1)

    C = 0.5
    list_C = []
    list_acc = []
    max_acc = 0
    max_C = 0
    while C<12.0:
        predictive_model = svm.SVC(C=C, kernel='rbf', gamma='auto')
        predictive_model.fit(train_X, train_y)
        acc = cross_val_score(predictive_model, test_X, test_y, cv = 10)
        list_C.append(C)
        list_acc.append(acc)

    plt.plot(list_C, list_acc)
    plt.xlabel("Parameter C")
    plt.ylabel("acc")
    plt.title("Optimizing Parameter C")
    plt.show()

normal_training()
