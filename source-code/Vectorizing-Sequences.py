import re



def convert_txt_to_pure_seq(string_seqs):
    new_string = string_seqs
    new_string.strip("BBa")
    pattern = r"[a-z]+"
    pure_seqs = re.findall(pattern, string_seqs)
    pure_seq = list(set(pure_seqs))
    return pure_seq




def positive_seq_list():
    f_positive_seq = open('Positive_seqs_with_ref_(101).txt', 'r')
    seq_positive = f_positive_seq.read()
    list_seq_positive = convert_txt_to_pure_seq(seq_positive)
    f_positive_seq.close()
    return list_seq_positive

def constitutive_seq_list():
    f_constitutive_seq = open('Constitutive_seqs_with_ref_(49).txt', 'r')
    constitutive_seq = f_constitutive_seq.read()
    list_seq_constitutive = convert_txt_to_pure_seq(constitutive_seq)
    list_seq_constitutive.remove('a')
    f_constitutive_seq.close()
    return list_seq_constitutive

def Negative_seq_list():
    f_Negative_seq = open('Negative_seqs_with_ref(130).txt', 'r')
    Negative_seq = f_Negative_seq.read()
    list_seq_negative = convert_txt_to_pure_seq(Negative_seq)
    list_seq_negative.remove('a')
    f_Negative_seq.close()
    return list_seq_negative

# 刚想到，可以用matplotlib查看一下序列长度的大致分布情况，然后考虑下一步

"""
def max_len(list_seq_positive, list_seq_constitutive, list_seq_negative):
    max_positive, max_constitutive, max_negative = 0, 0, 0
    for i in range(len(list_seq_positive)):
        max_positive = max(max_positive, len(list_seq_positive[i]))
    for i in range(len(list_seq_constitutive)):
        max_constitutive = max(max_constitutive, len(list_seq_constitutive[i]))
    for i in range(len(list_seq_negative)):
        max_negative = max(max_negative, len(list_seq_negative[i]))
    max_len = max(max_constitutive, max_negative, max_positive)
    print max_len
"""
