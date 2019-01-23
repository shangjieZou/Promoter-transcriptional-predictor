import re
import matplotlib.pyplot as plt


# This block is used to transform .txt into lists with pure sequences
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




# This block is for screening seqs which are oversize and plot distribution of seq length
def plot_distrbution_of_seq_len(list_seq_positive, list_seq_constitutive, list_seq_negative):
    def accumulate_seq_len(list_seq):
        x = []
        for i in range(len(list_seq)):
            x.append(len(list_seq[i]))
        return x
    x_positive = accumulate_seq_len(list_seq_positive)
    x_constitutive = accumulate_seq_len(list_seq_constitutive)
    x_negative = accumulate_seq_len(list_seq_negative)
    plt.hist(x_positive,histtype='stepfilled', alpha=0.4, density=6, bins=25, color="navy", label="Positive-regulated promoter")
    plt.hist(x_constitutive, histtype='stepfilled', alpha=0.4, density=6, bins=25, color="green", label="Constitutive promoter")
    plt.hist(x_negative, histtype='stepfilled', alpha=0.4, density=6, bins=25, color="darkgoldenrod", label="Repressible promoter")
    plt.ylabel("Frequency/Class Interval")
    plt.xlabel("Length of Sequences")
    plt.title("Distribution of Sequences' Length")
    plt.legend()
    plt.show()

def litering_oversize_seqs(list_seq):
    litered_seq_list = []
    for i in range(len(list_seq)):
        if len(list_seq[i]) <= 200:
            # if liter len_seq<=200, len(positive)=57, len(constitutive)=42, len(negative)=102
            litered_seq_list.append(list_seq[i])
    return litered_seq_list

def main_plot_distrbution():
    list_seq_positive = positive_seq_list()
    list_seq_constitutive = constitutive_seq_list()
    list_seq_negative = Negative_seq_list()
    liter_list_seq_positive = litering_oversize_seqs(list_seq_positive)
    liter_list_seq_constitutive = litering_oversize_seqs(list_seq_constitutive)
    liter_list_seq_negative = litering_oversize_seqs(list_seq_negative)
    plot_distrbution_of_seq_len(liter_list_seq_positive, liter_list_seq_constitutive, liter_list_seq_negative)


