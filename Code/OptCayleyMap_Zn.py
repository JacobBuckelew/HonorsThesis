""" 

OptCayleyMap_Zn.py is a large program that encompasses three different search algorithms to find an optimal
Cayley map embedding for a given Z_n group. The user will input the size of the group and be given the option
to choose a brute force search, hill climbing search, or stochastic hill climbing search

Jacob Buckelew
Fall 2021
Rollins College
"""

import math, random
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


# unique_perm generates the next largest unique permutation of a list of numbers
# returns a list of integers
def unique_perm(rho):
    # i will be used to iterate through the list of numbers starting from the end to find a rho[i] such that rho[i] > rho[i - 1]

    i = len(rho) - 1

    # Now loop through list starting from the end
    # 2 conditions must satisfy: there rho[i - 1] > rho[i] and that i > 1 since we dont care about comparing the 1 at the beginning
    # to the following number
    while((rho[i - 1] > rho[i] and i > 1)):
        i = i - 1
    
    # if i reached 0 after last iteration of the while loop, then the list is in descending order and is the greatest permutation possible
    # return false to tell program that the last unique permutation was already found
    if( i <= 1):
        return False
    

    # Now that we have found an index that satisfies rho[i] > rho[i - 1]
    # We need to search again from the end of the list for the largest rho[j] such that rho[j] > rho[i - 1]
    j = len(rho) - 1

    # We'll loop on one conditions: rho[j] <= rho[i - 1]
    # there must exist some rho[j] in this part of the subarray so we don't have to worry about j > 0
    while((rho[j] <= rho[i - 1])):
        j = j - 1
    
    # Now we will need to swap rho[i - 1] and rho[j]

    temp = rho[i - 1]
    rho[i - 1] = rho[j]
    rho[j] = temp


    # Finally we will need to reverse the sequence of elements rho[i] -> rho[length - 1]

    # Loop until i > j
    j = len(rho) - 1
    while(i < j):
        temp = rho[i]
        rho[i] = rho[j]
        rho[j] = temp
        i = i + 1
        j = j - 1

    
    return rho


def brute_search(n):
    # n - 1 elements in our rho
    elements = n - 1

    # keep track of the permutation that generates most optimal genus(based on Ringel and Young theorem)
    optimal_rho = list

    # also find the optimal genus according to Ringel and Young by calling optimal_genus(n)
    opt_genus = optimal_genus(n)


    # also store a string that will be used as input into a recursive function to find permutations of the
    # elements greater than 1 in the group of elements. 

    #We'll start with the rho in ascending order from 2...n-1

    rho = []
    for i in range(1,n):
        rho.append(i)


    # loop will be a stopping condition boolean

    loop = True
    genus = 0
    genus = calculate_genus(rho)
    # we'll assume that this genus is lowest for comparison reasons
    lowest_genus = genus
    low_rho = []
    low_rho = rho.copy()
    i = 1
    while((loop == True)):
        rho = unique_perm(rho)
        if(rho == False):
            loop = False
        # analyze rho here
        # keep track of a running lowest genus or find the optimal one according to ringel and young
        if(loop == True):
            genus = calculate_genus(rho)
            if(genus < lowest_genus):
                print("*Found lower genus*")
                print("rho:", rho)
                print("genus:", genus)
                lowest_genus = genus
                low_rho = rho.copy()
            i = i + 1
            if(i % 1000000 == 0):
                print("Program has checked " + str(i) + " permutations")
                print(rho)

    
    
    print("Solution found after " + str(i) + " iterations")

    return low_rho, lowest_genus



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
    
        

def stochastic_hc(current_permutation, genus, neighbors, opt_genus):

    best_permutation = current_permutation

    best_genus = genus

     # loop through the hill_climbing algorithm
    # keep running algorithm until we don't have any ways to get to another permutation
    # the algorithm will break from the loop and return the best rho and its genus
    # once the current permutation's neighbors do not have a lower genus


    # stochastic hill climbing requires looping for a definite interval of time
    k = 1
    while(k < 501):
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


    print("Best running genus after hill climber search was " + str(best_genus) + " with rho: " + str(best_permutation) + "\n")
    return best_permutation, best_genus


def standard_hc(current_permutation, genus, neighbors, opt_genus):

    starting_permutation = current_permutation
    # loop through the hill_climbing algorithm
    # keep running algorithm until we don't have any ways to get to another permutation
    # the algorithm will break from the loop and return the best rho and its genus
    # once the current permutation's neighbors do not have a lower genus
    
    previous_permutations = []



    while(len(neighbors) != 0):
        # keep track of next permutation to move to and its genus for comparison reasons
        next_permutation = None
        next_genus = 1000
        # iterate through all neighbors and compare each of their genuses to find 
        # most optimal neighbor
        for i in range(len(neighbors)):
            # for each neighbor check to see if they have been visited previously
            prev_perm = False
            for j in range(len(previous_permutations)):
                if(previous_permutations[j] == neighbors[i]):
                    prev_perm = True
            if(prev_perm != True):
                neighbor_genus = calculate_genus(neighbors[i])
                #print("Neighbor " + str(neighbors[i]) + " has genus " + str(neighbor_genus) )
                if(neighbor_genus <= next_genus):
                    next_permutation = neighbors[i]
                    next_genus = neighbor_genus
                    #if(next_genus == opt_genus):
                        #return next_permutation, next_genus
            else:
                continue
                #print("Neighbor #" + str(i) + " is a previously visited permutation")
        
        # then compare the neighboring genus to the current permutations genus
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
        #if(current_permutation == starting_permutation):
            #print("HI")
            #return current_permutation, genus
        neighbors = find_neighbors(current_permutation)


# hc_search will perform a local search using the hill climbing algorithm
# to find an optimal cayley map genus for the given Zn group

# will return a list representing rho and the genus of rho as an int
def hc_search(group, response):

    n = group
    opt_genus = optimal_genus(n)

    # before beginning the search, we will randomly generate a rho list
    # of length n - 1

    rho = generate_random_rho(n)

    # now begin the hill climbing algorithm using the randomly generated rho
    # define a current_permutation to keep track of which permutation
    # we are currently processing and comparing each neighboring permutation with

    current_permutation = rho

    # keep track of where we start
    starting_permutation = current_permutation


    # first evaluate the current permutation

    genus = calculate_genus(current_permutation)

    # If this rho has the optimal genus, then return this rho and its genus
    if(genus == opt_genus):
        print("*Starting genus was already optimal*\n")
        return current_permutation, genus

    neighbors = find_neighbors(current_permutation)
    
    # do a normal hill climb
    if(response == 2):
        return standard_hc(current_permutation, genus, neighbors, opt_genus)
    elif(response == 3):
        return stochastic_hc(current_permutation, genus, neighbors, opt_genus)

def perform_search(n, response):

    # conduct a brute force search
    if(response == 1):
        optimal_rho, genus = brute_search(n)
        print("The most optimal rho found by the brute force search was " + str(optimal_rho) + " which has a genus of " + str(genus) + "\n")


    elif(response == 2 or response == 3):
        # call hc_search(n) which will carry the bulk of the Hill climbing algorithm
        # returns the best rho the algorithm could find and its genus

        # 5 hill climbers will be simulated 

        best_hc_genus = 1000
        best_hc_rho = []

        for i in range(1,6):
            print("Beginning Hill Climber #" + str(i) + "'s search\n")
            best_rho, rho_genus = hc_search(n, response)

            if(rho_genus < best_hc_genus):
                best_hc_genus = rho_genus
                best_hc_rho = best_rho
            print("Hill Climber #" + str(i) + " found an optimal rho of " + str(best_rho) + " with genus " + str(rho_genus) + "\n")

        print("Most optimal rho found by the hill climber method: " + str(best_hc_rho) + "\n")
        print("Rho's genus: " + str(best_hc_genus) + "\n")
    
    else:
        print("Error in user input. Type a 1, 2, or 3 to indicate which method to use.\n")
    


# MAIN
continue_search = True


while(continue_search == True):
    # Ask for user input to decide on which Zn group we are working with
    print("Input an n to indicate which Z_n group you would like to calculate an optimal genus for:")
    n = int(input())

    print("Ringel and Young optimal genus for Z_" + str(n) + " is " + str(optimal_genus(n)) + "\n")


    # Ask user for what kind of search method they would like to use:
    # Brute force, Hill climbing or Stochastic Hill climb
    print("Input 1 for a brute force search\n")
    print("Input 2 for a hill climbing search\n")
    print("Input 3 for a stochastic hill climbing search\n")

    response = int(input())

    perform_search(n, response)

    print("Would you like to continue? Type 1 to do a new search or 2 to stop the program")
    
    continue_response = int(input())

    # end program
    if(continue_response == 2):
        continue_search = False
