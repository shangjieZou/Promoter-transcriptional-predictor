from Try_new_vectoring_method_with_dict import main_convert
import numpy as np
from sklearn.model_selection import cross_val_score,train_test_split
from sklearn import svm
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score

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
#y_list_positive = 89*[1]
#y_list_negative = 203 * [2]

#y_list = y_list_positive+y_list_constitutive+y_list_negative
y_list = y_list_constitutive + y_list_inducible
y = np.array(y_list).reshape(375)


#def normal_training():
train_X,test_X,train_y,test_y = train_test_split(X,y,test_size=0.2,random_state=1)

C = 0.5
list_C = []
list_acc = []
while C<48:
    #promoter_svm = svm.SVC(kernel='rbf', gamma="auto", C=C)
    #promoter_svm.fit(train_X, train_y)
    #predict_y = promoter_svm.predict(test_X)
    #acc = accuracy_score(test_y, predict_y)
    predictive_model = svm.SVC(C=C, kernel='rbf', gamma='auto')
    predictive_model.fit(train_X, train_y)
    acc_cross = cross_val_score(predictive_model, test_X, test_y, cv = 7)
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
plt.ylabel("ACC")
plt.title("Optimizing Parameter C")
plt.show()

