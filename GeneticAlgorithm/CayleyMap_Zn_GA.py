"""
CayleyMap_Zn_GA.py finds an optimal cayley map embedding for a given Z_n group using a genetic algorithm.


Jacob Buckelew
Fall 2021
Rollins College

"""

import math

# Define a population class to encapsulate all of the generations in a population
# as well as best_genus for each generation, average_genus for each generation, as
# well as the most optimal running genus and its permutation throughout all generations
class Population:
    def __init__(self, members):
        self.members = members

    
    


class Individual:
    def __init__(self, permutation):
        self.permutation = permutation
        print(self.permutation)
    
    def set_genus(self, genus):
        self.genus = genus

    
    def calculate_genus(self, rho):
        print(rho)
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
        

a = Individual([1,3,2])
a.set_genus(a.calculate_genus(a.permutation))


print(a.genus)


    




