from bwt_structures import *
from compsci260lib import *

def find(query, bwt_data):
    """
    Given a query sequence and a series of data structures
    containing various information about the reference genome,
    return a list containing all the locations of the query
    sequence in the reference genome. 
    """
    
    bwt, suffix_array, ranks, counts = bwt_data
      
    length = len(bwt)
    results = []
    

#     print "length: ", length
#     print "bwt: ", bwt
#     print "sorted bwt: ", sorted(bwt)
#     print "suffix_array: ", suffix_array
#     print "ranks: ", ranks
#     print "counts: ", counts

    start_range_F = 0
    realign = 0
    query = query[::-1]
    for i in range(len(query)):
#         print query[i]
        F_index = counts.get(query[i]) + 1
        start_range_F = F_index + realign
        if i == 0: 
            num_letters = ranks.get(query[i])[-1]
        end_range_F = start_range_F + num_letters - 1 
        if i < len(query) - 1:
            next_letter = query[i+1]
#         print "start_range_F:", start_range_F
#         print "start_range_L:", start_range_F
        
        try:
            start_range_L = ranks.get(next_letter)[start_range_F]
            end_range_L = ranks.get(next_letter)[end_range_F]
        except:
            return []
        
        num_letters = end_range_L - start_range_L
        if bwt[start_range_F] == next_letter: 
            num_letters += 1
           
        x = 1 
        
        for j in range(start_range_F, end_range_F + 1):
            if bwt[j] == next_letter and x == 1:
                realign = ranks.get(next_letter)[j] - 1 
                x = 0
            
#     print "start: ", start_range_F
#     print "end: ", end_range_F
        
    for i in range(start_range_F, end_range_F+1):
        results.append(suffix_array[i])
#     print sorted(results)
    return sorted(results)
    
if __name__ == '__main__':
    # example query sequence
    query_sequence = "AAACGA"
    # example reference sequence
    sequence = "AAAAAAAAACGATAGAGA"
    find(query_sequence, make_all(sequence))
    sequence = "AAAAAAAAACGATAGAGAAAAAAAAAACGATAGAGA"
    find(query_sequence, make_all(sequence))
    sequence = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    find(query_sequence, make_all(sequence))
    
    
    
