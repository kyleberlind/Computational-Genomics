from compsci260lib import *
import sys, random, os
import textwrap
from sys import exit
import random
from numpy.random import choice
import tidy_text
from itertools import count

# seq_length is the length of the sequence we will generate
def generate_HMM(hmm_file,seq_length):
    if not os.path.exists(hmm_file):
        print("Can't open HMM parameter file: %s" % hmm_file)
        return -1

    f = open(hmm_file, "r")

    # read the state names
    states = f.readline().strip().split()

    # read the initial probabilities
    initial_probs = f.readline().strip().split()

    # read the transition matrix
    transitions = [None for _ in range(len(states))]
    for i in range(0, len(states)):
        matrix_row = f.readline().strip().split()
        transitions[i] = matrix_row
        
    # read the input alphabet
    input_alphabet = f.readline().strip().split()
    
    # read the input matrix
    inputs = [None for _ in range(len(states))]
    for i in range(0, len(states)):
        matrix_row = f.readline().strip().split()
        inputs[i] = matrix_row

    f.close()

#     print "matrix row", matrix_row
#     print "inputs: ", inputs
#     print "transitions: ", transitions
#     print "states: ", states
#     print "initial_probs: ", initial_probs
    
    
    solve_HMM(1000, states, transitions, initial_probs, inputs)
    
    
def solve_HMM(num_outtext, states, transitions, initial_probs, inputs):
    curr_state = []
    HMM_genome = []
    previous_state = choice(states, 1, p = initial_probs)[0]
    for i in range(0,1000):
        state = get_state(previous_state, states, transitions)
        previous_state = state
        if state == "state1":
            curr_state.append("1")
        if state == "state2":
            curr_state.append("2")
        HMM_genome.append(generate_next_nucleotide(state, inputs))    
     
    HMM_genome = "".join(HMM_genome)
    state_list = "".join(curr_state)
    
    print state_list
    print HMM_genome
    
    print "\n"
    print "State\t\tNumber of positions spent in that state"
    
    j = 0 
    i = 0 

    length_num_positions =0
    length_num_positions1 = 0 
    length_num_positions2 = 0 
    count2 = 0
    count1 = 0
    while i < len(state_list):
        s = state_list[i]
        num_positions = 0
        while j < len(state_list) and state_list[j] == s: 
            num_positions += 1
            j+=1
        print "%s\t\t%s" % (s, num_positions)
        if s == "1": 
            count1 += num_positions
            length_num_positions1 += 1
        if s == "2": 
            count2 += num_positions
            length_num_positions2 += 1
        i = j
        num_positions = 0
        length_num_positions += 1
        

    avg_state_length1 = count1/float(length_num_positions1)
    avg_state_length2 = count2/float(length_num_positions2)
    print "average length in state1: ", avg_state_length1
    print "average length in state2: ", avg_state_length2

def generate_next_nucleotide(state, inputs):
    '''get a random next character following seed from 
    the text dictionary'''
    inputs2 = [0.2, 0.3, 0.3, 0.2]
    inputs1 = [0.3, 0.2, 0.2, 0.3]
    if state == "state1": 
        next_char = choice(inputs[0], 1, p = inputs1)[0] #here I use numpy to weight the choices based on a probability distribution

    if state == "state2":
#         inputs[1] = ['.2', '.3', '.3', '.2']
        next_char = choice(inputs[0], 1, p = inputs2)[0]# change this
#         print copyinputs
    return next_char

def get_state(previous_state, states, transitions):
    
    state1_trans_prob = transitions[0]
    state2_trans_prob = transitions[1]
    
#     print" state1_trans_prob: ", state1_trans_prob
#     print" state2_trans_prob: ", state2_trans_prob
#     print states
    if previous_state == "state1":
        next_state = choice(states, 1, p= state1_trans_prob)[0]
#         print choice(states, 1, p = state1_trans_prob)[0]
#         print state1_trans_prob
    if previous_state == "state2": 
        next_state = choice(states, 1, p= state2_trans_prob)[0]
    
    return next_state
    
if __name__ == '__main__':
    # you can change the parameters if necessary
    generate_HMM("../HMM.txt",1000)
