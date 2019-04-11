import sys
from math import log, exp
from compsci260lib import get_fasta_dict

def posterior_decoding(f_in_file, f_hmm_file):
    # read the state names and number
    states = f_hmm_file.readline().split()
    K = len(states)
    
    # read the initial probabilities
    probs = f_hmm_file.readline().split()
    initial_probs = [float(prob) for prob in probs]
    
    # read the transition matrix
    transitions = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [float(trans_prob) for trans_prob in matrix_row_arry]
        transitions[i] = matrix_row
    
    # read the emission symbols
    emission_symbols = f_hmm_file.readline().split()
    
    # read the emission probability matrix
    emit_probs = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [float(emit_prob) for emit_prob in matrix_row_arry]
        emit_probs[i] = matrix_row
    
    f_hmm_file.close()
    
    seq_dict = get_fasta_dict(input_file)
    emit_str = seq_dict.values()[0]  # there's only 1
    
    print "Done reading sequence of length ", len(emit_str)
    
    # Run the forward algorithm
    forward = run_forward(states, initial_probs, transitions, emission_symbols, emit_probs, emit_str)
    
    # Run the backward algorithm
    backward = run_backward(states, initial_probs, transitions, emission_symbols, emit_probs, emit_str)
    
    # Calculate the posterior probabilities
    # Initializing the posterior 2D matrices
    posterior = [[float(0) for _ in range(K)] for _ in range(len(emit_str))]
    for i in range(len(emit_str)):
        # Did not normalize the probabilities (i.e., did not divide by P(X)),
        # because we will only use these probabilities to compare
        # posterior[i][0] versus posterior[i][1]
        for k in range(K):
            posterior[i][k] = forward[i][k] + backward[i][k]

    # Print out the decoded results
    
#     for post in posterior:
#         max_p
#         for k in range(K):
    trace = []
    for post in posterior:
    
        max_post= float("-inf")
        max_state = 0
        for k in range(K):
            if post[k] > max_post:
                max_post = post[k]
                max_state = k
        trace.append(max_state)
#     print trace[1627030:1627050]
#     print trace [1627041-1:1627042+2]

#     display_trix(trace)
    catch_top_ten(trace)
    catch_bottom_ten(trace)
        

def run_forward(states, initial_probs, transitions, emission_symbols, emit_probs, emit_str):
    """Calculates the forward probability matrix"""

    K = len(states)
    L = len(emit_str)
    
#     print "emission probs", emit_probs
#     print "transitions", transitions
    
    forward = [[float(0) for _ in range(K)] for _ in range(L)]

    # Initialize
    emit_index = get_emit_index(emit_str[0], emission_symbols)
    for k in range(K):
        forward[0][k] = log(initial_probs[k]) + log(emit_probs[k][emit_index])

    # Iterate
#     print len(forward)
    for i in range(1, L):
        emit_index = get_emit_index(emit_str[i].upper(), emission_symbols)
        for k in range(K):
            my_values = []
            for j in range(K):
                ajk = transitions[j][k]
                ej = emit_probs[k][emit_index]
                fk = forward[i-1][j]
                vali = log(ajk) + log(ej) + fk
                my_values.append(vali)
            forward[i][k] = sum_vals(my_values[0], my_values[1])
      
        # Compute the forward probabilities for the states
    return forward        

def sum_vals(val1, val2):
    ret = val1 + log(1 + exp(val2 - val1))
    return ret

def run_backward(states, initial_probs, transitions, 
    emission_symbols, emit_probs, emit_str):
    """Calculates the backward probability matrix"""

    K = len(states)
    L = len(emit_str)

    backward = [[float(0) for _ in range(K)] for _ in range(L)]

    # Initialize
    for k in range(K):
        backward[L-1][k] = log(1)  # which is zero, but just to be explicit...
    
    # Iterate
    for i in range(L-2, -1, -1):
        emit_index = get_emit_index(emit_str[i+1].upper(), emission_symbols)
        for j in range(K):
            my_values = []
            for k in range(K):
                ajk = transitions[j][k]
                ej = emit_probs[k][emit_index]
                fk = backward[i+1][k]
                vali = log(ajk) + log(ej) + fk
                my_values.append(vali)
            backward[i][j] = sum_vals(my_values[0], my_values[1])
    
        # Compute the backward probabilities for the states
    return backward   
     
def display_trix(traceback):
    i = 1
    j = 0
    count = 0
    print "segment\tstart\tstop\tstate\tlength in state "
    while (i<len(traceback)): 
        cur_state = traceback[i-1]
        while (j < len(traceback) and traceback[j] == cur_state):
            j+=1
        ind_start = i
        ind_stop = j 
        len_in_state = (j - i+1)
        if cur_state == 0: 
            print_cur_state = "state1"
        if cur_state == 1: 
            print_cur_state = "state2"
        print "%d\t%s\t%d\t%s\t%d" % (count, ind_start, ind_stop, print_cur_state, len_in_state)
        count += 1
        j+=1
        i=j    
    print "len trace back", count
    
def catch_top_ten(traceback):
    i = 1
    j = 0
    count = 0
    print "First Ten Segments"
    print "segment\tstart\tstop\tstate"
    while (i<len(traceback)): 
        cur_state = traceback[i-1]
        while (j < len(traceback) and traceback[j] == cur_state):
            j+=1
        ind_start = i 
        ind_stop = j 
        count += 1
        if cur_state == 0: 
            print_cur_state = "state1"
        if cur_state == 1: 
            print_cur_state = "state2"
        if count <= 10: 
            print "%d.\t%s\t%d\t%s" % (count, ind_start, ind_stop, print_cur_state)
    
        j+=1
        i=j 
        
def catch_bottom_ten(traceback):
    i = 1
    j = 0
    count = 0
    tcount = 0 
    print "Last Ten Segments"
    print "segment\tstart\tstop\tstate"
    while (i<len(traceback)): 
        cur_state = traceback[i-1]
        while (j < len(traceback) and traceback[j] == cur_state):
            j+=1
            
        tcount += 1
        j+=1
        i=j 
    i = 1
    j = 0
    state2count = 0
    while (i<len(traceback)): 
        cur_state = traceback[i-1]
        if cur_state+1 == 2:
            state2count += 1
        while (j < len(traceback) and traceback[j] == cur_state):
            j+=1
        ind_start = i 
        ind_stop = j 
        count += 1
        if cur_state == 0: 
            print_cur_state = "state1"
        if cur_state == 1: 
            print_cur_state = "state2"
        if count > tcount - 10: 
            print "%d.\t%s\t%d\t%s" % (count, ind_start, ind_stop, print_cur_state)
        j+=1
        i=j 
        
    print "Structural RNA Regions: ", state2count

def get_emit_index(input_val, alphabet):
    """does a linear search to find the index of a character"""
    for i in range(len(alphabet)):
        if alphabet[i] == input_val:
            return i
    sys.exit("Could not find character " + input_val)


if __name__ == '__main__':
    hmm_file = "../HMM.methanococcus.txt"
    input_file = "../bacterial.genome.fasta"

    f_in_file = open(input_file)
    f_hmm_file = open(hmm_file)
    
    if f_in_file is None:
        sys.exit("Can't open HMM file: " + hmm_file)
    if f_hmm_file is None:
        sys.exit("Can't open file: " + input_file)
    
    posterior_decoding(f_in_file, f_hmm_file)
