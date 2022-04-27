""" 

CayleyMap_Zn_HC.py applies a hill climbing approach to solving the problem of minimizing 
of a Cayley Map embedding for the finite cyclic group Zn 

Jacob Buckelew
Fall 2021
Rollins College
"""

import math, random, timeit
import numpy as np
from numpy.random import rand


# optimal_genus calculates the optimal genus for a complete graph Kn embedding
# Based on the Ringel and Young Theorem
# returns an int
def optimal_genus(n):
    genus = math.ceil(((n-3)*(n-4))/12)
    return genus



# generate_rho(n) will take an integer as an argument to represent which Zn group
# we are making a rho for
# will return a list that represents a randomly generated permutation for rho


def generate_random_rho(n):

    # create a list to represent rho
    # this rho will always begin with 1

    random_rho = [1]

    # Use numpy to generate a random permutation of numbers from 2 to n -1
    random_rho = random_rho + list(np.random.permutation(np.arange(2,n)))


    return random_rho


# calculate_genus(rho) will return an integer representing the genus of a given
# rho, which is inputted in the form of a list
def calculate_genus(rho):
    # calculating a genus  involves first finding the lambdas

    # First make a dictionary to hold all of the lambda(key) values given some key which represents a number in the rho
    # 1...n-1
    # This map holds the lambda(x) values at their corresponding keys
    lambda_map = {}
    n = len(rho)

     # generate a copy of the lambda_map dictionary(compare_map) which will be used to keep track of which
    # elements of rho have already been used
    # each element will map to 1s and then change these values to 0 as we come across the element while
    # figuring out our lambdas

    compare_map = {}

    # edge case to consider when filling the lambda(x) dictionary)
    # be careful of an element that exists at the end of a rho

    last_element = rho[len(rho) - 1]

    # find its inverse 

    last_ele_inv = (n + 1) - last_element


    for i in range(n):

        # cover edge case
        if(rho[i] == last_element):
            lambda_map[last_ele_inv] = rho[0]
        else:
        # fill in lambda normally
            lambda_map[(n + 1) - rho[i]] = rho[i + 1]

        # fill in the compare_map key, value pair

        compare_map[rho[i]] = 1
    

    # We now want to loop from 1 up to n + 1
    # and iteratively find the faces using the lambda_map
    # edges and vertices can be calculated using simple formulas

    # We know that edges
    edges = ((n * (n + 1))/(2))

    vertices = n + 1

    faces = 0

    euler_char = 0

    genus = 0

    comparison = [0] * n
    values = list(compare_map.values()) 
    j = 1

    for i in range(1, n + 1):
        # j will be used to keep track of which value we move to next
        # will be
        j = i
        counter = 0
        m = j
        face_lambda = 0
        mult = 0
        values = list(compare_map.values())
        if(values == comparison):
            break


        if(compare_map[j] == 1):

            while(compare_map[j] != 0):
                compare_map[j] = 0
                counter = counter + 1
                if(compare_map[lambda_map[j]] != 0):
                    m = m + lambda_map[j]
                else:
                    break
                j = lambda_map[j]
        else:
            continue


        # now that we have m and n we can do some math to find face_lambda
        # and then iterate the faces variable to the number of faces generated from the lambda
        # we just created

        # Now use a formula derived from the relationship between gcd(m,n), m, n, and lcm(m,n)
        # to find the multiplicity of the lambda that was found

        mult = ((n + 1)/(math.gcd(m, n + 1)))
        face_lambda = mult * counter

        faces = faces + (((n+1) * (counter))/(face_lambda))
        
    # now calculate the euler characteristic
    euler_char = vertices - edges + faces 


    # now calculate the genus

    genus = int((euler_char - 2)/(-2))


    return genus


# find_neighbors(rho) will return a list of all permutations 
# that are potential next permutations given a rho which is in the form of a list



def find_neighbors(rho):

    # In order to find neighbors, we will do some swaps between adjacent elements of the permutation
    
    i = 0

    neighbors = []

    while(i < len(rho)):

        # edge case when we reach the end of rho, 
        # we will want to swap the first with the last
        if(i == len(rho) - 1):
            temp = rho[i]
            rho[i] = rho[0]
            rho[0] = temp
            neighbors.append(rho.copy())
            temp = rho[i]
            rho[i] = rho[0]
            rho[0] = temp

        # else swap normally between adjacent elements
        else:
            temp = rho[i]
            rho[i] = rho[i + 1]
            rho[i + 1] = temp
            neighbors.append(rho.copy())
            temp = rho[i]
            rho[i] = rho[i + 1]
            rho[i + 1] = temp

        i = i + 1

    
    return neighbors



def random_transition(neighbors, neighbor_genuses, n):

    # create a list to store all of the candidates
    # each element will be an altered genus from neighbor_genuses
    # this will scale each genus so that some genuses are more likely
    # to be chosen then others to reach the optimal genus in the long run
    candidates = {}

    # keep track of each interval to figure out which candidate gets picked
    # by the random number generation process
    intervals = []

    # keep track of the size of the search space so we can generate a random number from
    # (0, search_size)
    search_size = 0

    previous_interval = 0


    # for i in range(len(neighbor_genuses)):
    #     # subtract the genus from a number that scales with group size
    #     # and then square the genuses
    #     genus = (math.ceil(((n)** 2)/4) - neighbor_genuses[i]) ** 2
    #     # print("New genus: ", genus)
    #     if(genus > 0):
    #         # append interval
    #         interval = previous_interval + genus
    #         intervals.append(interval)
    #         previous_interval = interval
    #         search_size = search_size + genus
    #         # save information at this candidates location in the dictionary
    #         candidates[interval] = [neighbor_genuses[i], neighbors[i]]
    #     else:
    #         continue
    
    # print("Candidates: ", candidates)
    # print('Search size:', search_size)
    # print(intervals)

    genus_mappings = {}
    collections = {}
    genus_intervalmap = {}
    genuses = []

    interval = 0

    for i in range(len(neighbor_genuses)):
        old_genus = neighbor_genuses[i]
        genus = (math.ceil(((n)** 2)/4) - neighbor_genuses[i]) ** 2

        if(genus > 0):

            if(genus not in genus_mappings.keys()):
                genus_mappings[genus] = old_genus
                genuses.append(genus)
                interval = previous_interval + genus
                intervals.append(interval)
                previous_interval = interval
                candidates[genus] = []
                candidates[genus].append([neighbor_genuses[i],neighbors[i]])
                collections[interval] = []
                collections[interval].append(neighbors[i])
                genus_intervalmap[genus] = interval
                search_size = search_size + genus
            else:
                candidates[genus].append([neighbor_genuses[i], neighbors[i]])
                collections[genus_intervalmap[genus]].append(neighbors[i])
            
        else:
            continue


    # for i in range(len(neighbor_genuses)):
    #     genus = (math.ceil(((n)** 2)/4) - neighbor_genuses[i]) ** 2
        
    #     if ((genus > 0)):
    #         print("interval:", interval)
    #         if(genus not in genuses):
    #             genuses.append(genus)
    #         if():
    #             interval = previous_interval + genus
    #             print("new interval: ", interval)
    #             genus_mappings[genus] = interval
    #             print(genus_mappings)
    #             candidates[interval] = []
    #             intervals.append(interval)
    #             previous_interval = interval
    #             search_size = search_size + genus
    #         candidates[interval].append([neighbor_genuses[i], neighbors[i]])
    #     else:
    #         continue


    # Now randomly generate a number in the search spaces

    random_int = random.randint(1, search_size)


    previous_interval = 0
    for i in range(len(intervals)):
        # found interval
        if((random_int > previous_interval) and (random_int <= intervals[i])):
            #print(random_genus)

            # randomly choose a list from the list of values at the intervals[i] key

            interval_values = candidates.get(genuses[i])

            transition_neighbor = random.choice(interval_values)


            return transition_neighbor[1], transition_neighbor[0]
        previous_interval = intervals[i]
    
        





# hc_search will perform a local search using the hill climbing algorithm
# to find an optimal cayley map genus for the given Zn group

# will return a list representing rho and the genus of rho as an int
def hc_search(group, solution_start):

    n = group
    opt_genus = optimal_genus(n)



    # before beginning the search, we will randomly generate a rho list
    # of length n - 1

    rho = generate_random_rho(n)

    print("Beginning search at rho:", rho)

    # now begin the hill climbing algorithm using the randomly generated rho
    # define a current_permutation to keep track of which permutation
    # we are currently processing and comparing each neighboring permutation with

    current_permutation = rho

    # keep track of best running genus and rho as we go
    # sometimes we will move away from the best running genus so we must keep track of it

    best_permutation = current_permutation

    # first evaluate the current permutation

    genus = calculate_genus(current_permutation)
    
    best_genus = genus

    print("Starting Rho has genus ", genus)

    # If this rho has the optimal genus, then return this rho and its genus
    if(genus == opt_genus):
        print("*Starting genus was already optimal*")
        solution_end = timeit.default_timer()
        solution_time = solution_end - solution_start
        print("Solution found in " + str(solution_time) + " seconds")
        return current_permutation, genus

    neighbors = find_neighbors(current_permutation)
    


    # loop through the hill_climbing algorithm
    # keep running algorithm until we don't have any ways to get to another permutation
    # the algorithm will break from the loop and return the best rho and its genus
    # once the current permutation's neighbors do not have a lower genus


    # stochastic hill climbing requires looping for a definite interval of time
    k = 1
    while(k < 500):
        next_permutation = None
        transition = False
        next_genus = 1000
        # randomly choose an order to  to examine using the random sample function

        size_neighbors = len(neighbors)

        neighbor_genuses = []


        for i in range(size_neighbors):
            neighbor_genus = calculate_genus(neighbors[i])
            #print("Neighbor " + str(neighbors[i]) + " has genus " + str(neighbor_genus) )
                
             # Decide whether to accept transitioning to this neighbor or not

            neighbor_genuses.append(neighbor_genus)


        next_permutation, next_genus = random_transition(neighbors, neighbor_genuses, n)
            
        #print("Next Permuation:", next_permutation)
        #print("Next genus: ", next_genus)

        # then compare the neighboring genus to the current permutations genus
        # if none of the neighboring genuses are lower, then just return where we are at
        # also return if we have traveled in a circle

        #print("Genus: ", next_genus)
        #print("Best_genus: ", best_genus)
        if(next_genus < best_genus):
            best_genus = next_genus
            #print("new best genus: ", best_genus)
            best_permutation = next_permutation
            #print("new best permutation: ", best_permutation)

        if(best_genus == opt_genus):
            #print("Found optimal genus " + str(best_genus) + " after " + str(k) + " iterations")
            return next_permutation, best_genus

        # move to next permutation and loop again with new neighbors
        # save the previous permutation to ensure we don't loop endlessly
        current_permutation = next_permutation
        #print("Current permutation: ", current_permutation)
        #print("New current permutation: ", current_permutation)
        neighbors = find_neighbors(current_permutation)

        k = k + 1


    print("Best running genus after search was " + str(best_genus) + " with rho: " + str(best_permutation))
    return best_permutation, best_genus

# MAIN

# Ask for user input to decide on which Zn group we are working with
print("Input an n to indicate which Z_n group you would like to calculate an optimal genus for:")
n = int(input())

print("Ringel and Young optimal genus for Z_" + str(n) + " is " + str(optimal_genus(n)) + "\n")


# call hc_search(n) which will carry the bulk of the Hill climbing algorithm
# returns the best rho the algorithm could find and its genus

# 5 hill climbers will be simulated 

best_hc_genus = 1000
best_hc_rho = []
for i in range(1,6):
    solution_start = timeit.default_timer()
    print("Beginning Hill Climber #" + str(i) + "'s search")
    best_rho, rho_genus = hc_search(n, solution_start)

    if(rho_genus < best_hc_genus):
        best_hc_genus = rho_genus
        best_hc_rho = best_rho
    
    solution_end = timeit.default_timer()
    solution_time = solution_end - solution_start
    print(rho_genus)
    print("Solution found in " + str(solution_time) + " seconds ")
    print("Hill Climber #" + str(i) + " found an optimal rho of " + str(best_rho) + " with genus " + str(rho_genus) + "\n")

print("Best rho: ", best_hc_rho)
print("rho genus: ", best_hc_genus)