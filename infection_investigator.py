from bwt_structures import *
from read_aligner import *
from compsci260lib import *

bacteria_proportion_count_dict = {"Bacteroides_ovatus": 0,
                        "Bacteroides_thetaiotaomicron": 0,
                        "Bifidobacterium_longum": 0,
                        "Eubacterium_rectale": 0,
                        "Lactobacillus_acidophilus": 0,
                        "Peptoniphilus_timonensis": 0,
                        "Prevotella_copri": 0,
                        "Roseburia_intestinalis": 0,
                        "Ruminococcus_bromii": 0,
                        "Vibrio_cholerae": 0}

def reverse_complement(seq):
    """
    Returns the reverse complement of the input string.
    """
    comp_bases = {'A': 'T',
                  'C': 'G',
                  'G': 'C',
                  'T': 'A'} 
    rev_seq = list(seq)
    rev_seq = rev_seq[::-1]
    rev_seq = [comp_bases[base] for base in rev_seq] 
    return ''.join(rev_seq)
    
def align_patient_reads():
    
    """ this function build several dictionaries to store values and indecies associtated with comparing reads
    to the genomes of different bacteria. bacteria_dict contains the genomes (both the forward and reverse strands)
    mapped to the bateria. bacteria data dict contains the makeall() data for each bacterias genome. bacteria count
    dict stored the number of times that a read matches a genome, and bacteria indecies dict stores the indicies that
    my find function generates. This function will print out the prevalences of the bacteria in each of the patients
    by calling calc proportions"""
    
    
    #NOTE: THERE IS A BUG IN THIS FUNCTION THAT PREVENTS ME FROM FINDING THE CORRET PREVALANCE OF CHOLERA IN ONE OF 
    #THE PATIENTS. 
    
    #PLEASE TAKE A LOOK AT IT, I THINK IT MIGHT BE SOMETHING SMALL. THINK ITS BECAUSE I USE NUMBER OF MATCHES IN THE GENOME,
    #NOT JUST THE NUMBER OF READS THAT MATCH A BACTERIA TO CALCULATE THE PREVELANCE

    patient_list = []
    
    patient1_dict = get_fasta_dict("patients/patient1.fasta")
    patient2_dict = get_fasta_dict("patients/patient2.fasta")
    patient3_dict = get_fasta_dict("patients/patient3.fasta")
    patient_list.append(patient1_dict)
    patient_list.append(patient2_dict)
    patient_list.append(patient3_dict)
    
    bacteria_dict = {"Bacteroides_ovatus": [get_fasta_dict("reference_genomes/Bacteroides_ovatus.fasta").values()[0], reverse_complement(get_fasta_dict("reference_genomes/Bacteroides_ovatus.fasta").values()[0])],
                    "Bacteroides_thetaiotaomicron": [get_fasta_dict("reference_genomes/Bacteroides_thetaiotaomicron.fasta").values()[0], reverse_complement(get_fasta_dict("reference_genomes/Bacteroides_thetaiotaomicron.fasta").values()[0])],
                    "Bifidobacterium_longum": [get_fasta_dict("reference_genomes/Bifidobacterium_longum.fasta").values()[0], reverse_complement(get_fasta_dict("reference_genomes/Bifidobacterium_longum.fasta").values()[0])],
                    "Eubacterium_rectale": [get_fasta_dict("reference_genomes/Eubacterium_rectale.fasta").values()[0], reverse_complement(get_fasta_dict("reference_genomes/Eubacterium_rectale.fasta").values()[0])],
                    "Lactobacillus_acidophilus":[get_fasta_dict("reference_genomes/Lactobacillus_acidophilus.fasta").values()[0], reverse_complement(get_fasta_dict("reference_genomes/Lactobacillus_acidophilus.fasta").values()[0])],
                    "Peptoniphilus_timonensis":[get_fasta_dict("reference_genomes/Peptoniphilus_timonensis.fasta").values()[0], reverse_complement(get_fasta_dict("reference_genomes/Peptoniphilus_timonensis.fasta").values()[0])],
                    "Prevotella_copri": [get_fasta_dict("reference_genomes/Prevotella_copri.fasta").values()[0], reverse_complement(get_fasta_dict("reference_genomes/Prevotella_copri.fasta").values()[0])],
                    "Roseburia_intestinalis":[get_fasta_dict("reference_genomes/Roseburia_intestinalis.fasta").values()[0], reverse_complement(get_fasta_dict("reference_genomes/Roseburia_intestinalis.fasta").values()[0])],
                    "Ruminococcus_bromii":[get_fasta_dict("reference_genomes/Ruminococcus_bromii.fasta").values()[0], reverse_complement(get_fasta_dict("reference_genomes/Ruminococcus_bromii.fasta").values()[0])],
                    "Vibrio_cholerae": [get_fasta_dict("reference_genomes/Vibrio_cholerae.fasta").values()[0], reverse_complement(get_fasta_dict("reference_genomes/Vibrio_cholerae.fasta").values()[0])]}
    
    bacteria_data_dict = {}
    
    for bacteria, sequences in bacteria_dict.items():
        bacteria_data_dict[bacteria] = [make_all(sequences[0]), make_all(sequences[1])]
    

    
    bacteria_count_dict = {"Bacteroides_ovatus": 0,
                           "Bacteroides_thetaiotaomicron": 0,
                           "Bifidobacterium_longum": 0,
                           "Eubacterium_rectale": 0,
                           "Lactobacillus_acidophilus": 0,
                           "Peptoniphilus_timonensis": 0,
                           "Prevotella_copri": 0,
                           "Roseburia_intestinalis": 0,
                           "Ruminococcus_bromii": 0,
                           "Vibrio_cholerae": 0}
    

    for i in range(len(patient_list)):
        for read_num, read in patient_list[i].items():
            num_non_empty = 0
            bacteria_indecies_dict = {}
            for bacteria, bacteria_seq in bacteria_dict.items():
                if bacteria not in bacteria_indecies_dict.keys():
                    bacteria_indecies_dict[bacteria] = []
                find_run_forward = find(read, bacteria_data_dict[bacteria][0])
                if find_run_forward != []:
                    bacteria_indecies_dict[bacteria].extend(find_run_forward) 
                    num_non_empty += 1
                find_run_backward = find(read, bacteria_data_dict[bacteria][1])
                if find_run_backward != []:
                    bacteria_indecies_dict[bacteria].extend(find_run_backward)
                    num_non_empty += 1
            if num_non_empty > 1:
                continue
            if  num_non_empty == 1:
                for k, v in bacteria_indecies_dict.items():
                    if v != []:
                        bacteria_count_dict[k] += len(v)# counts the length of the list at that indicie and adds it to the count

        print "Patient" + str(i + 1) + " Prevalance:\n"
        calc_proportions(bacteria_count_dict)
        print "\n"
    


def calc_proportions(bacteria_count_dict):
    """calculates the bacterial prevalences in a patient"""
    total_count = 0
    for bacteria, count in bacteria_count_dict.items():
        total_count += count
    for bacteria, count in bacteria_count_dict.items():
        bacteria_proportion_count_dict[bacteria] = float(count)/total_count
    for bacteria, proportion in bacteria_proportion_count_dict.items():
        print bacteria, proportion

if __name__ == '__main__':
    align_patient_reads()