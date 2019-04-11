
import sys
from math import log
from compsci260lib import *

def viterbi_decoding(input_file, f_hmm_file):
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
    
    #take the logs of the transition probabilities

    
    # read the emitted symbols
    emitted_symbols = f_hmm_file.readline().split()

    # read the emission probability matrix
    emit_probs = [None for _ in range(K)]
    for i in range(K):
        matrix_row_arry = f_hmm_file.readline().split()
        matrix_row = [float(emit_prob) for emit_prob in matrix_row_arry]
        emit_probs[i] = matrix_row

    #take the logs of the emission probabilities
            
    f_hmm_file.close()
    
    seq_dict = get_fasta_dict(input_file)
    emit_str = seq_dict.values()[0]  #there's only 1

    print "Done reading sequence of length ", len(emit_str)
    
    # Create Viterbi table and traceback table
    viterbi = [[0 for _ in range(len(emit_str))] for _ in range(K)]   
    pointers = [[0 for _ in range(len(emit_str))] for _ in range(K)] 

    # Initialize the first column of the matrix

    for i in range(K):
        in_index = get_emit_index(emit_str[0].upper(), emitted_symbols)
        viterbi[i][0] = log(emit_probs[i][in_index]) + log(initial_probs[i])
        
    for probability in emit_probs:
        for i in range(len(probability)):
            probability[i] = log(probability[i])
            
    
    for probability in transitions: 
        for i in range(len(probability)):
            probability[i] = log(probability[i])
    
    # Build the matrix column by column
    for j in range(1, len(emit_str)):#zero index is already filled by this point
        in_index = get_emit_index(emit_str[j].upper(), emitted_symbols)
        for i in range(K):
            el = emit_probs[i][in_index]
#             akl = transition[i][
#             vl = max([(viterbi[x][j-1] + transitions[x][i]) for x in range(K)])
            
            max_vl = float('-inf')
            state_vl = 0
            for x in range(K):
                vl = viterbi[x][j-1] + transitions[x][i]
                if vl > max_vl: 
                    max_vl = vl
                    state_vl = x + 1
               
            viterbi[i][j] = (el + max_vl)
            pointers[i][j] = state_vl
    

            # Compute the entries viterbi[i][j] and pointers[i][j]
            # NEG_INF is float('-inf')

       
    max_v = float('-inf')
    for i in range(K): 
        pointers[i] = pointers[i][1:]
        if viterbi[i][-1] > max_v: 
            max_v = viterbi[i][-1]
            ind_max_v = i
   
    traceback = []
    i = len(pointers[0])
    state2 = ind_max_v + 1
    while i > 0:
        state1 = state2
        traceback.append(state1)
        i -= 1
        state2 = pointers[state1-1][i]

#reverse the list
    traceback = traceback[::-1]
#     display_trix(traceback)
    catch_top_ten(traceback)
    catch_bottom_ten(traceback)
    # Traceback code goes here:
    
def display_trix(traceback):
    i = 0
    j = 0

    print "start\tstop\tstate\tlength in state "
    while (i<len(traceback)): 
        cur_state = traceback[i]
        while (j < len(traceback) and traceback[j] == cur_state):
            j+=1
            
        ind_start = i + 1
        ind_stop = j + 1
        
        len_in_state = (j - i)
        print "%s\t%d\t%s\t%d" % (ind_start, ind_stop, cur_state, len_in_state)
        j+=1
        i=j 

         
#     print "start\tstop\tstate "
#     for  in :
#         print "%s\t%d\t%s" % (key, txt_dict[key][0], txt_dict[key][1] )

def catch_top_ten(traceback):
    i = 0
    j = 0
    count = 0
    print "First Ten Segments"
    print "segment\tstart\tstop\tstate"
    while (i<len(traceback)): 
        cur_state = traceback[i]
        while (j < len(traceback) and traceback[j] == cur_state):
            j+=1
        ind_start = i + 1
        ind_stop = j + 1
        count += 1
        if count <= 10: 
            print "%d.\t%s\t%d\t%s" % (count, ind_start, ind_stop, cur_state)
    
        j+=1
        i=j 
def catch_bottom_ten(traceback):
    i = 0
    j = 0
    count = 0
    tcount = 0 
    print "Last Ten Segments"
    print "segment\tstart\tstop\tstate"
    while (i<len(traceback)): 
        cur_state = traceback[i]
        while (j < len(traceback) and traceback[j] == cur_state):
            j+=1
        ind_start = i + 1
        ind_stop = j + 1
        tcount += 1
        j+=1
        i=j 
        
    i = 0
    j = 0
    state2count = 0
    while (i<len(traceback)): 
        cur_state = traceback[i]
        if cur_state == 2:
                state2count += 1
        while (j < len(traceback) and traceback[j] == cur_state):
            j+=1
        ind_start = i + 1
        ind_stop = j + 1
        count += 1
        if count > tcount - 10: 
            print "%d.\t%s\t%d\t%s" % (count, ind_start, ind_stop, cur_state)
        j+=1
        i=j 
    print "Structural RNA Regions: ", state2count
    
def get_emit_index(input_val, alphabet):
    for i in range(len(alphabet)):
        if alphabet[i] == input_val:
            return i
    
    sys.stderr("Could not find character " + input_val)

if __name__ == '__main__':
    hmm_file = "../HMM.methanococcus.txt"
    input_file = "../bacterial.genome.fasta"
    
    f_in_file = open(input_file)
    f_hmm_file = open(hmm_file)
    
    if f_in_file is None:
        sys.exit("Can't open HMM file: " + hmm_file)
    if f_hmm_file is None:
        sys.exit("Can't open file: " + input_file)
    
    viterbi_decoding(input_file, f_hmm_file)
