from compsci260lib import *
import textwrap
from sys import exit
import random
from numpy.random import choice
import tidy_text


def solve_language(file_name,order,num_outtext):
    # Input of the function is the order of the Markov Model
    # text file processed by tidy.text.py
    # e.g. file_name = 'tidy.heart.of.darkness.txt'

    # Read the contents of the file
    f = open(file_name, 'r')

    if f is None:
        print "Can't open " + file_name
    else:
        contents = f.read()
        f.close()
        contents = contents.replace("\n", "")
        contents = contents.replace("\r", "")

    # This dictionary will store all the data needed to estimate the Markov model:
    txt_dict = {}

    # Step 1: Count up occurrences of the various k-tuples in the training text
    print "Building dict of k-tuples and their counts..."
    build_dict(contents, order, txt_dict)

    # Step 2: Collect the counts necessary to estimate transition probabilities
    print "Collecting counts to estimate transition probabilities..."
    collect_counts(contents, order, txt_dict)
    display_dict(txt_dict)

    # Step 3: Generate artificial text from the trained model
    for _ in range(num_outtext):
        seed = contents[0:order]
        M = 1000
        next_character = ""
        text = seed

        print "\nOne version of the story is:"

        # Display my story
        for _ in range(M):
            next_character = generate_next_character(seed, txt_dict)
            text += next_character
            seed = seed[1:] + next_character

        text_list = textwrap.wrap(text, 72)
        text = "\n".join(text_list)
        print text

def display_dict(txt_dict):
    print "key\tcount\tfollowers "
    for key in txt_dict:
        print "%s\t%d\t%s" % (key, txt_dict[key][0], txt_dict[key][1] )

def build_dict(contents, k, txt_dict):
    '''builds a dictionary of k ordered tuples, where tuples are matched 
    to the number of times they occur in a training text.'''
    
    print contents
    
    for i in range(len(contents) -k):# might break
        if not contents[i:i+k] in txt_dict: 
            txt_dict[contents[i:i+k]] = 0
        txt_dict[contents[i:i+k]] += 1 
    return txt_dict

#update(k,v)
def collect_counts(contents, k, txt_dict):
    '''collects the counts of the letters following each
     key in text dictionary built in build_dict'''
    
    for i in range(len(contents) -k):
        count = txt_dict[contents[i:i+k]]
        nuc = contents[i:i+k]
        if type(txt_dict[contents[i:i+k]]) == int:
            txt_dict[nuc] = (count, {})
        if contents[i+k] not in txt_dict[nuc][1]:             
            txt_dict[nuc][1][contents[i+k]] = 0
        txt_dict[nuc][1][contents[i+k]] += 1 
    
    print txt_dict
    return txt_dict
        
        
def generate_next_character(seed, txt_dict):
    '''get a random next character following seed from 
    the text dictionary'''
    
    possibility_dict = txt_dict[seed][1]
    options = [key for key in possibility_dict.keys()]
    total = float(sum([value for value in possibility_dict.values()]))
    weightings = [value/total for value in possibility_dict.values()]#turn occurances into probabilities by dividing by total occurances
#     print "options: ", options
#     print "weightings: ", weightings
    next_char = choice(options, 1, p = weightings)[0] #here I use numpy to weight the choices based on a probability distribution

#     print choice(options, 1, weightings)
    return next_char
    

if __name__ == '__main__':

    file_name = "../tidy.heart.of.darkness.txt"
    order = 1
    num_outtext = 1
    solve_language(file_name,order,num_outtext)
