""""
CayleyMap_Zn determines an optimal cayley map genus for the finite cyclic group Z_n


Jacob Buckelew
Fall 2021
Rollins College
"""



import math
import matplotlib.pyplot as plt


# optimal_genus calculates the optimal genus for a complete graph Kn embedding
# Based on the Ringel and Young Theorem
# returns an int
def optimal_genus(n):
    genus = math.ceil(((n-3)*(n-4))/12)
    return genus


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

# calculate_genus solves for the lambdas of a given rho and then calculates the genus for that rho
# returns an int representing the genus, as well as the permutation

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

    print(rho)
    print(genus)

    return genus

# graph_results will produce a plot that shows the overall distribution
# of genus across the rhos in a certain groyp
def graph_results(xt, y, x):

    fig = plt.figure(figsize=(10,10), dpi=75 )
    plt.xticks(x, xt, rotation="vertical")
    plt.plot(y)
    plt.show()



# find_genus is the main function controlling the analysis of each rho permutation
# It will call unique_perm(rho) to find the next rho and also calculate lambdas to find a genus for each rho
# returns an optimal permutation in the form of a list
def find_genus(n):
    # n - 1 elements in our rho
    elements = n - 1

    # keep track of the permutation that generates most optimal genus(based on Ringel and Young theorem)
    optimal_rho = list

    # also find the optimal genus according to Ringel and Young by calling optimal_genus(n)
    opt_genus = optimal_genus(n)
    print("The Ringel and Young optimal genus for Z_" + str(n) + " is ", opt_genus)


    # also store a string that will be used as input into a recursive function to find permutations of the
    # elements greater than 1 in the group of elements. 

    #We'll start with the rho in ascending order from 2...n-1

    rho = []
    for i in range(1,n):
        rho.append(i)

    
    # xt will save the rhos for the x axis of the genus graph
    xt = []

    # y will save the y values correpsponding to a rho's genus
    y = []

    # loop will be a stopping condition boolean

    loop = True
    genus = 0
    genus = calculate_genus(rho)
    # we'll assume that this genus is lowest for comparison reasons
    lowest_genus = genus
    low_rho = []
    low_rho = rho.copy()
    x_label = rho.copy()


    xt.append(x_label)
    y.append(genus)

    i = 1
    while((loop == True)):
        rho = unique_perm(rho)
        if(rho == False):
            loop = False
        # analyze rho here
        # keep track of a running lowest genus or find the optimal one according to ringel and young
        if(loop == True):
            genus = calculate_genus(rho)
            x_label = rho.copy()
            xt.append(x_label)
            y.append(genus)
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

    x = list(range(0,i))
    print(xt)
    print(x)
    graph_results(xt, y, x)

    
    
    print("Solution found after " + str(i) + " iterations")

    return low_rho, lowest_genus
       

# MAIN

print("Input an n to indicate which Z_n group you would like to calculate an optimal genus for:")
n = int(input())
optimal_rho, genus = find_genus(n)

print("The most optimal rho found by the program was " + str(optimal_rho) + " which has a genus of " + str(genus))