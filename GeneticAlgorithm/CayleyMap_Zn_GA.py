"""
CayleyMap_Zn_GA.py finds an optimal cayley map embedding for a given Z_n group using a genetic algorithm.


Jacob Buckelew
Fall 2021
Rollins College

"""

import math, random

# Object Classes

# Define a population class to encapsulate all of the generations in a population
# as well as best_genus for each generation, average_genus for each generation, as
# well as the most optimal running genus and its permutation throughout all generations
class Population:

    opt_genus = 90
    def __init__(self, members):
        self.members = members
    
    def set_members(self, members):
        self.members = members

    def set_ringelgenus(self, ringel_genus):
        self.ringelgenus = ringel_genus

    def set_best_genus(self, best_genus):
        self.best_genus = best_genus
    
    def set_avg_genus(self, avg_genus):
        self.avg_genus = avg_genus
    
    
    def set_best_rho(self, rho):
        self.best_rho = rho

    
    # metrics will iterate through the population to find avg_genus, best_genus, opt_genus, and opt_rho
    # for each generation
    def metrics(self):

        ringel_genus = optimal_genus(n)
        # keep track of average, best, most optimal running rho and its genus
        # population_members will store individuals

        avg_genus = 0
        best_genus = 1000
        best_rho = []



        
        # iterate through the members
        for member in self.members:
            # find individual's genus

            genus = member.calculate_genus(member.permutation)

            member.set_fitness(genus)

            print(member.fitness)

            # set best_genus if member fitness is less than running best genus
            if(member.fitness < best_genus):
                # update population's optimal rho and genus
                best_rho = member.permutation
                print("new best genus: ", self.best_genus)
                opt_genus = genus
                self.best_genus = genus
            
            # Print if we see a rho that satisfies ringel definition, this will be most optimal genus and rho
            if(member.fitness == self.ringelgenus):
                print("Ringel and Young optimal genus found")
                print("Genus: ", member.fitness)
                print("Rho: ", member.permutation)
                self.best_genus = genus

            # add up genus to running total to find the average genus in the generation

            avg_genus = avg_genus + genus
            


        print("Done")


    
    

# Define an individual class
class Individual:
    def __init__(self, permutation):
        self.permutation = permutation
        #print(self.permutation)
    
    def set_fitness(self, genus):
        self.fitness = genus

    def get_fitness(self):
        return self.fitness

    
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
        


# Functions

# optimal_genus calculates the optimal genus for a complete graph Kn embedding
# Based on the Ringel and Young Theorem
# returns an int
def optimal_genus(n):
    genus = math.ceil(((n-3)*(n-4))/12)
    return genus


# define random_permutation(n) to generate random permutations for a given group size(n)
# returns a permutation in the form of a list
def random_permutation(n):
    permutation = []
    
    # fill permutation with 1, 2, ... n-1
    i = 1
    while i < n:
        permutation.append(i)
        i = i + 1

    # use random.shuffle to obtain a random permutation
    random.shuffle(permutation)
    return permutation

    



def ga_search(n):

    population_size = 100
    # population_members will store individuals

    # keep track of average, best, most optimal running rho and its genus
    # population_members will store individuals
    population_members = []


    for i in range(population_size):
        # initialize individuals
        member = Individual(random_permutation(n))

        # add individual to the population of individuals
        population_members.append(member)
    
    population = Population(population_members)
    
    population.set_members(population_members)
    population.set_best_genus(1000)


    ringel_genus = optimal_genus(n)
    population.set_ringelgenus(ringel_genus)

    # calculate metrics of the initial population

    population.metrics()



   

# MAIN

print("Input an n: ")
n = int(input())

ga_search(n)










    




