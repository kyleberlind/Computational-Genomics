import sys, random
from compsci260lib import *
from math import exp


#2 a.
G = 3*10**6
R = 4*10**4
L = 450

avg_emp_coverage = []
avg_num_uncovered_nucleotides = [] 
avg_num_contigs = []
avg_avg_length_contigs = [] 

def simulate():
    """takes in the length of the genome G, the number of reads R, and the length of each read L and
    Computes the coverage C which is expected number of times each nucleotide in the genome has been 
    sequenced during the procedure, expected number of nucleotides in the genome that remain unsequenced 
    during this procedure, an expression for the expected number of nucleotides in the genome that remain 
    unsequenced during this procedure, the expected number of contigs reported by such an algorithm, and
    the expected length of each contig """
    
    genome = []
    
#  print end_read
    for i in range(G):
        genome.append(0)
    for j in range(R):
        start_read = random.randint(0,G)
        end_read = start_read + L
        if end_read >= len(genome):
                end_read = len(genome)
        for k in range(start_read, end_read):
            genome[k] += 1
   
    zero_count = 0 
    for l in range(G):
        if genome[l] == 0:
            zero_count += 1
            
        
   
    emp_coverage = float(sum(genome))/len(genome) #this will yeild 5.99999 because we are cutting off the edge cases, making a slightly smaller coverage than there should be
    print "empirical coverage: ", emp_coverage
    print "number of uncovered nucleotide", zero_count

    count_contigs(genome)
    print "\n"
    
    avg_emp_coverage.append(emp_coverage)
    avg_num_uncovered_nucleotides.append(zero_count)

    
    
    
def count_contigs(genome):
    """ counts the number of contigs in a given sequence by 
    findin breaks in the sequence and appending the contiguous
    chuncks to an array. returns the avg num and avg len 
    of the contigs it finds."""
    contigs = []
    j = 0
    for i in range(len(genome)):
        if genome[i] == 0:
            if j != i:
                contigs.append(i-j)
            j = i + 1
        
    num_contigs = len(contigs)
    avg_len_contigs = sum(contigs)/num_contigs
    

#     contigs = genome
#     contigs = ''.join(str(i) for i in contigs)
#     contigs = contigs.split("0")
#     contigs = [contig for contig in contigs if contig != ""]
#     print contigs[:20]
#     num_contigs = len(contigs)
#     avg_len_contigs = 0
#     for contig in contigs:
#         avg_len_contigs += len(contig)
#     avg_len_contigs = avg_len_contigs/num_contigs
#     
    print "number of contigs: ", num_contigs
    print "avgerage length of a contig: ", avg_len_contigs
            
    avg_num_contigs.append(num_contigs)
    avg_avg_length_contigs.append(avg_len_contigs)

    
def run_simulation():
    """ runs the simulation code to print out each
    run in a readable format"""
    
    
    for i in range(20):
        print "Run" + str(i) + ": "
        simulate()
    
if __name__ == '__main__':
    run_simulation()
    
    avg_emp_coverage = sum(avg_emp_coverage)/len(avg_emp_coverage)
    avg_num_uncovered_nucleotides = sum(avg_num_uncovered_nucleotides)/len(avg_num_uncovered_nucleotides)
    avg_num_contigs = sum(avg_num_contigs)/len(avg_num_contigs)
    avg_avg_length_contigs =sum(avg_avg_length_contigs)/len(avg_avg_length_contigs)
    
    print "average empirical coverage: ", avg_emp_coverage
    print "average number of uncovered nucleotides: ", avg_num_uncovered_nucleotides
    print "average number of of contigs: ", avg_num_contigs
    print "average of average length of contigs: ", avg_avg_length_contigs
    
    