import itertools



def eumerate():
    """enumerates the possible permutations of given states"""
    
    outcomes = set(itertools.combinations_with_replacement("FLFL", 4))
    
#     print len(outcomes)
    list = []
    for thing in outcomes: 
        list.append("".join(thing))
    list.append("LFLF")
    list = sorted(list)
    
    count = 0
    probability = 0
    probabilities = []
    for element in list:
        for i in range(len(element)):
            if i == 0:
                probabilities.append("(0.5)")
                probability = 0.5
            else:
                if element[i] == "F" and element[i-1] == "F":
                    probabilities.append("(0.95)")
                    probability *= 0.95
                elif element[i] == "F" and element[i-1] == "L":
                    probabilities.append("(0.1)")
                    probability *= 0.1
                elif element[i] == "L" and element[i-1] == "L":
                    probabilities.append("(0.9)")
                    probability *= 0.9
                else:
                    probabilities.append("(0.05)")
                    probability *= 0.05
        count += 1
        print str(count) + ".", element + ":", "P(" + element + ")", "=", "".join(probabilities), "=", probability
        probabilities = []
        probability = 0 

if __name__ == '__main__':
    eumerate()