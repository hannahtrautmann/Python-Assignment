import pandas as pd
import sys
import plotnine as p9

def get_sequence(filename):
    """
    Summary line: Reads in a file

    Extended description: Takes an input filename and outputs a string

    Parameters:
    filename : sequence filename as a string

    Return:
    string : entire sequence as a string
    """

    with open(filename, 'r') as current_file:
        text = current_file.read()
    return(text)

def observed_kmer_dict_single_value(sequence, k):
    """
    Summary line: Generates a dictionary of kmers from a sequence of a given k

    Extended description:Takes a sequence, determines possible ks that could be made, cycles through those adding each unique kmer to a dictionary.
    Adds a count to the dictionary if not a new kmer.

    Parameters:
    sequence : a string of nucleotides
    k (int) : a positive integer used to determine size of kmers

    Return:
    dictionary : a dictionary that contains all unique kmers in the given sequence, with counts of how many found in that sequence
    """
    assert k >0, "Please specify k > 0"
    assert k <=len(sequence), "Please specify k less than or equal to length of sequence"
    observed_kmer_dict = {}
    total_kmers = len(sequence) - k + 1
    for i in range(total_kmers):
        kmer=sequence[i:i+k]
        if kmer not in observed_kmer_dict:
            observed_kmer_dict[kmer] = 1
        else:
            observed_kmer_dict[kmer] += 1
    return observed_kmer_dict

def observed_kmers_list(sequence):
    """
    Summary line: Generates a list of observed kmers for a given sequence

    Extended description:Takes a sequence and loops through reading kmers,
    adding a value to the dictionary if it is a new.
    Then counts values in each dictionary to determine unique kmers for each k.

    Parameters:
    sequence : a string of nucleotides

    Return:
    list : a list that contains number of all observed kmers for each possible k
    """
    observed = []
    for obk in range(1,len(sequence)+1):
        observed_kmer_dict = {}
        total_kmers = len(sequence) - obk + 1
        for i in range(total_kmers):
            kmer=sequence[i:i+obk]
            if kmer not in observed_kmer_dict:
                observed_kmer_dict[kmer] = 1
            else:
                observed_kmer_dict[kmer] += 1
        observed.append(len(observed_kmer_dict))
    return observed

def possible_kmers_list(sequence):
    """
    Summary line: Generates a list of possible kmers for a given sequence

    Extended description:Takes a sequence and determines if 4^k is smaller than the length of the sequence -k + 1 for any given k.
    Whichever value is smaller for each k, that is the possible kmers. Adds this value to a list.

    Parameters:
    sequence : a string of nucleotides

    Return:
    list : a list that contains number of possible kmers for each possible k
    """
    possible = []
    for i in range(1,len(sequence)+1):
        if (len(sequence) - i + 1) > (4**i):
            possible.append(4**i)
        else:
            possible.append(len(sequence) - i + 1)
    return possible

def kmers_df(sequence):
    """
    Summary line: Generates a datraframe of observed and possible kmers for a given sequence

    Extended description:Takes a sequence and calculates possible k values based on sequence length.
    Creates a dictionary of all observed kmers.
    Then counts number of entries in that dictionary to determine unique observerd kmers.
    Next calculates possible kmers for any sequence of that length.
    Creates a dictionary of the observed, possible, and k values.
    Turns dictionary into a dataframe using pandas

    Parameters:
    sequence : a string of nucleotides

    Return:
    dataframe : a dateframe that contains list of k values, observed kmers, and possible kmers
    """
    observed = []
    for obk in range(1,len(sequence)+1):
        observed_kmer_dict = {}
        total_kmers = len(sequence) - obk + 1
        for i in range(total_kmers):
            kmer=sequence[i:i+obk]
            if kmer not in observed_kmer_dict:
                observed_kmer_dict[kmer] = 1
            else:
                observed_kmer_dict[kmer] += 1
        observed.append(len(observed_kmer_dict))
    possible = []
    for i in range(1,len(sequence)+1):
        if (len(sequence) - i + 1) > (4**i):
            possible.append(4**i)
        else:
            possible.append(len(sequence) - i + 1)
    k = list(range(1,len(sequence)+1))
    obs_poss_dict = {'k':k, 'observed_kmers':observed, 'possible_kmers':possible}
    import pandas as pd
    kmers_df = pd.DataFrame(obs_poss_dict, columns = ['k','observed_kmers', 'possible_kmers'])
    return kmers_df


def make_plots(kmers_df,filename):
    """
    Summary line: Creates plot of kmer proportion_kmers

    Extended description: Takes dataframe and generates a scatterplot with k on the x-axis and proportion of observed over possible kmers on the y axis
    Then saves this plot based ont he filename

    Parameters:
    dataframe : dataframe that contains counts of possible ks, possible kmers, and observed kmers_df for a sequence
    filename: name of file that has source sequence in it

    Return:
    plot: plot saved to pdf with name of sequence file and description
    """
    proportion_kmers = p9.ggplot(data=kmers_df, mapping = p9.aes(x= 'k', y = 'observed_kmers/possible_kmers')) + p9.geom_point()
    proportion_kmers.save(filename+"_proportion_kmers.pdf")

def linguistic_complexity(sequence):
    """
    Summary line: Calculates linguistic complexity for a given sequence

    Extended description:Takes a sequence and first generates a dictionary of observed kmers.
    Then counts number of entries to determine unique observerd kmers.
    Next calculates possible kmers for any sequence of that lenght.
    Creates a dictionary of the observed, possible, and k values.
    Turns dictionary into a dataframe using pandas
    Adds all observed observed kmers and adds all possible kmers
    Calculates linguistic complexity by dividing total possible by total observed

    Parameters:
    sequence : a string of nucleotides

    Return:
    decimal : a single value between 0 and 1 representing linguistic complexity
    """
    assert sequence != '', "Please specify sequence that is not empty"
    observed = []
    for obk in range(1,len(sequence)+1):
        observed_kmer_dict = {}
        total_kmers = len(sequence) - obk + 1
        for i in range(total_kmers):
            kmer=sequence[i:i+obk]
            if kmer not in observed_kmer_dict:
                observed_kmer_dict[kmer] = 1
            else:
                observed_kmer_dict[kmer] += 1
        observed.append(len(observed_kmer_dict))
    possible = []
    for i in range(1,len(sequence)+1):
        if (len(sequence) - i + 1) > (4**i):
            possible.append(4**i)
        else:
            possible.append(len(sequence) - i + 1)
    k = list(range(1,len(sequence)+1))
    obs_poss_dict = {'k':k, 'observed_kmers':observed, 'possible_kmers':possible}
    import pandas as pd
    kmers_df = pd.DataFrame(obs_poss_dict, columns = ['k','observed_kmers', 'possible_kmers'])
    total_observed = kmers_df['observed_kmers'].sum()
    total_possible = kmers_df['possible_kmers'].sum()
    linguistic_complexity = total_observed/total_possible
    return linguistic_complexity

#this "main" section is run only if the entire script is being run on command line. It is not used for pytesting.
if __name__ == '__main__':
    myfile = sys.argv[1] #put filename as first argument after name of this code on command line
    myseq = get_sequence(myfile) #generates a squence from the file
    my_kmers_df = kmers_df(myseq) #creates dataframe with k, possible and observed kmers
    my_linguistic_complexity = linguistic_complexity(myseq) #calculates linguistic complexity
    print(my_kmers_df.iloc[0:5]) #returns first five lines of dataframe so appears on command line
    print('The linguistic complexity is', my_linguistic_complexity) #returns linguistic complexity so appears on command line
    make_plots(my_kmers_df,myfile) #generates pdf file with plot of proportion kmers
