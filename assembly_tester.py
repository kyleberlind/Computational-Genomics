from compsci260lib import *
from operator import itemgetter
import operator
def assembly_tester(reads,contigs,super_contigs):
    super_contig_dict = get_fasta_dict(super_contigs)
    contig_dict = get_fasta_dict(contigs)
    read_dict = get_fasta_dict(reads)
    
    
    check1 = read_in_set(read_dict,contig_dict)
    check2 = check_distance(read_dict,contig_dict)
    check3 = ordered_contig_check(read_dict,contig_dict, super_contig_dict)
    
    if not check1: 
        print "check 1 failed"
    if not check2:
        print "check 2 failed"
    if not check3:
        print "check 3 failed"
    if check1 and check2 and check3:
        print "assembly is consistent with the reads"
        
    
# check that every read appears somewhere in the set of contigs in the correct orientation.
def read_in_set(read_dict,contig_dict):
    print "running first check:"
    reads_in_contigs_list = []
    reads_in_contigs = False
    for read, sequence in read_dict.items():
        for contig, ref_seq in contig_dict.items():
            if sequence in ref_seq: 
#                 print read
                reads_in_contigs = True
        reads_in_contigs_list.append(reads_in_contigs)
        if not reads_in_contigs: 
            print read, "not found!"
        reads_in_contigs = False
            
    overall_check = True
    for status in reads_in_contigs_list:
#         print status
        overall_check = overall_check and status
#     print "Check 1:", overall_check
    return overall_check
    
def check_distance(read_dict,contig_dict):
#check that whenever a mated pair of reads appears in the same contig, the distance 
#from the beginning of the first read to the end of the last read is 2000 + or - 20.
    
    contig_list = [value for value in contig_dict.values()]
    contains_list = [[] for i in contig_dict.keys()]
    #map the first contig in the dict to contig one, the second to contig two

    for read, sequence in read_dict.items():
        for i in range(len(contig_list)):
            if sequence in contig_list[i]:
                contains_list[i].append(read)
    pairs_list = []
    unmatched_contigs = []
    for item in contains_list:
        contig_index = contains_list.index(item)
        item.sort() 
        for i in range(len(item)):
            #check if the matched reads are in the same contig
            if item[i][:-1] + "a"  in item and item[i][:-1] + "b" in item and (item[i][:-1] + "a", item[i][:-1] + "b") not in pairs_list:
                pairs_list.append((item[i][:-1] + "a", item[i][:-1] + "b", contig_index))
            if item[i][:-1] + "a"  not in item and item[i][:-1] + "b" in item  or item[i][:-1] + "a" in item and item[i][:-1] + "b" not in item:
                unmatched_contigs.append((item[i][:-1] + "a", item[i][:-1] + "b"))
    #pairs list will have all of the reads that appear on the same contig, unmatched contigs will have all the reads that appear on different contigs
    print "\nrunning second check:"
    for pair in pairs_list:
        start = read_dict[pair[0]]
        end = read_dict[pair[1]]
        #get the referenced contig from the tuple
        ref_contig = contig_list[pair[2]]
        start_index = ref_contig.index(start) 
        end_index = ref_contig.index(end) + len(end)
        
        dist = end_index - start_index
        if dist<=1980 or dist >= 2020:
#             print "Check 2:", False
            return False
        
#     print "Check 2:", True
    return True

def ordered_contig_check(read_dict,contig_dict,super_contigs):

#check that whenever a mated pair of reads appears in the same contig, the distance 
#from the beginning of the first read to the end of the last read is 2000 + or - 20.
    
    contig_list = [value for value in contig_dict.values()]
    contains_list = [[] for i in contig_dict.keys()]
    #map the first contig in the dict to contig one, the second to contig two
    

    for read, sequence in read_dict.items():
        for i in range(len(contig_list)):
            if sequence in contig_list[i]:
                contains_list[i].append(read)
    unmatched_contigs = []
    for item in contains_list:
        item.sort() 
        for i in range(len(item)):
            #check if the matched reads are in different contig
            if item[i][:-1] + "a"  not in item and item[i][:-1] + "b" in item  or item[i][:-1] + "a" in item and item[i][:-1] + "b" not in item:
                unmatched_contigs.append((item[i][:-1] + "a", item[i][:-1] + "b"))
    
    print"\nrunning third check:\n"
    unmatched_contigs = set(unmatched_contigs)
    unmatched_contigs = list(unmatched_contigs)
    
    for pair in unmatched_contigs:
        ref_contig = super_contigs['super_contig0']
        start = read_dict[pair[0]]
        end = read_dict[pair[1]]
        start_index = ref_contig.index(start) 
        end_index = ref_contig.index(end) + len(end)
        
        if end_index - start_index < len(start) + len(end):
            print pair
            
#         print "start ", start_index, "end ", end_index
        
        if end_index < start_index or not (start in ref_contig and end in ref_contig):
#             print "Check 3: ", False
            return False
#     print "Check 3:", True
    return True
        
        
    
    
                
    
    
                
                
                
                
    
            
    
    

if __name__ == '__main__':
    assembly_tester("../paired.reads.fasta","../contigs_combo.fasta"  ,"../supercontig.fasta")
