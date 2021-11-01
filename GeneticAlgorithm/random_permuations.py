# This program simply exists for comparison purposes to gain an understanding of how well the genetic algorithm can
# actually perform against a random search.


import math, random

# define random_permutation(n) to generate random permutations for a given group size(n)
# returns a permutation in the form of a list
def random_permutation(n):
    permutation = []
    
    # fill permutation with 1, 2, ... n-1
    permutation.append(1)
    i = 2
    elements = []
    while i < n:
        elements.append(i)
        i = i + 1

    # use random.shuffle to obtain a random permutation
    random.shuffle(elements)

    permutation = permutation + elements

    # print(permutation)

    if(permutation == [1,6,15,7,17,16,5,18,2,9,13,14,11,4,10,12,8,3]):
        print("Found optimal *****************************")
    return permutation


# optimal_genus calculates the optimal genus for a complete graph Kn embedding
# Based on the Ringel and Young Theorem
# returns an int
def optimal_genus(n):
    genus = math.ceil(((n-3)*(n-4))/12)
    return genus


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


# MAIN

print("Input an n")

n = int(input())

best_genus = 100
best_rho = list

opt_genus = optimal_genus(n)


print("Ringel and Young optimal genus: ", opt_genus)

for i in range(100000):
    rho = random_permutation(n)
    genus = calculate_genus(rho)

    print(rho)
    print(genus)

    if(rho ==[1,6,15,7,17,16,5,18,2,9,13,14,11,4,10,12,8,3]):
        print("OPTIMAL********")
        break

    if(genus == best_genus):
        print(str(rho) + " also has genus " + str(best_genus))

    if(genus < best_genus):
        print(str(rho) + " has next lower genus of " + str(genus))
        best_rho = rho
        best_genus = genus
    
print("Best rho: " + str(best_rho) + " with genus of " + str(best_genus))




