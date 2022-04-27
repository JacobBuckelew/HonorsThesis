""" 

CayleyMap_Zn_HC.py applies a hill climbing approach to solving the problem of minimizing 
of a Cayley Map embedding for the finite cyclic group Zn 

Jacob Buckelew
Fall 2021
Rollins College
"""

import math, random, timeit
import numpy as np


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

    # keep track of where we start
    starting_permutation = current_permutation


    # first evaluate the current permutation

    genus = calculate_genus(current_permutation)

    print("Starting Rho has genus ", genus)

    # If this rho has the optimal genus, then return this rho and its genus
    if(genus == opt_genus):
        print("*Starting genus was already optimal*")
        solution_end = timeit.default_timer()
        solution_time = solution_end - solution_start
        print("Solution found in " + str(solution_time) + " seconds ")
        return current_permutation, genus

    neighbors = find_neighbors(current_permutation)
    


    # loop through the hill_climbing algorithm
    # keep running algorithm until we don't have any ways to get to another permutation
    # the algorithm will break from the loop and return the best rho and its genus
    # once the current permutation's neighbors do not have a lower genus
    
    previous_permutations = []



    while(len(neighbors) != 0):
        # keep track of next permutation to move to and its genus for comparison reasons
        next_permutation = None
        #print("All previous permutations:", previous_permutations)
        next_genus = 1000
        #print("Rho's Neighbors: ", neighbors)
        # iterate through all neighbors and compare each of their genuses to find 
        # most optimal neighbor
        for i in range(len(neighbors)):
            #print("Next_genus: ", next_genus)
            # for each neighbor check to see if they have been visited previously
            prev_perm = False
            for j in range(len(previous_permutations)):
                if(previous_permutations[j] == neighbors[i]):
                    #print("Found previous permutation:", neighbors[i])
                    prev_perm = True
            if(prev_perm != True):
                neighbor_genus = calculate_genus(neighbors[i])
                #print("Neighbor " + str(neighbors[i]) + " has genus " + str(neighbor_genus) )
                if(neighbor_genus <= next_genus):
                    next_permutation = neighbors[i]
                    #print("Next permutation", next_permutation)
                    next_genus = neighbor_genus
                    #if(next_genus == opt_genus):
                        #return next_permutation, next_genus
            else:
                print("Neighbor #" + str(i) + " is a previously visited permutation")
        
        # then compare the neighboring genus to the current permutations genus
        #print("Next lower genus:", next_genus)
        # if none of the neighboring genuses are lower, then just return where we are at
        # also return if we have traveled in a circle
        if((next_genus > genus)):
            return current_permutation, genus
        else:
            genus = next_genus

        # move to next permutation and loop again with new neighbors
        # save the previous permutation to ensure we don't loop endlessly
        previous_permutations.append(current_permutation)
        current_permutation = next_permutation
        #print("Current permutation: ", current_permutation)
        #print("Starting permutation:, ", starting_permutation)
        #if(current_permutation == starting_permutation):
            #print("HI")
            #return current_permutation, genus
        #print("New current permutation: ", current_permutation)
        neighbors = find_neighbors(current_permutation)


    


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
    print("Hill Climber #" + str(i) + " found an optimal rho of " + str(best_rho) + " with genus " + str(rho_genus) + "\n")
    solution_end = timeit.default_timer()
    solution_time = solution_end - solution_start
    print(rho_genus)
    print("Solution found in " + str(solution_time) + " seconds ")
print("Best rho: ", best_hc_rho)
print("rho genus: ", best_hc_genus)