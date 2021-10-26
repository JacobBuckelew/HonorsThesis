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
    best_genus = []
    avg_genus = []
    best_rho = []
    sum_fitness = 0
    members = list
    ringelgenus = 0


    def __init__(self, members):
        self.members = members
    
    def set_members(self, members):
        self.members = members

    def set_ringelgenus(self, ringel_genus):
        self.ringelgenus = ringel_genus

    def set_best_genus(self, best_genus):
        self.best_genus.append(best_genus)
    
    def set_best_rho(self, best_rho):
        self.best_genus.append(best_rho)
    
    def set_avg_genus(self, avg_genus):
        self.avg_genus.append(avg_genus)

    def set_sum_fitness(self, sum_fitness):
        self.sum_fitness = sum_fitness
    
    
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

            # If statement to determine best genus so far
            # set best_genus if member fitness is less than running best genus for this generation
            if(member.fitness < best_genus):
                # update population's optimal rho and genus
                best_rho = member.permutation
                best_genus = genus
            
            # Print if we see a rho that satisfies ringel definition, this will be most optimal genus and rho
            if(member.fitness == self.ringelgenus):
                #print("Ringel and Young optimal genus found")
                #print("Genus: ", member.fitness)
                #print("Rho: ", member.permutation)
                best_genus = genus
                best_rho = member.permutation

            

            # add up genus to running total to find the average genus in the generation

            avg_genus = avg_genus + genus
            
        sumfitness = avg_genus
        self.set_sum_fitness(sumfitness)
        avg_genus = int(avg_genus/len(self.members))
        # store best genus in list for this generation
        self.set_best_genus(best_genus)

        # store avg genus in list for this generation
        self.set_avg_genus(avg_genus)
        
        # store best rho in this generation
        self.set_best_rho(best_rho)


        #print("Best genus for this generation: ", best_genus)
        #print("best rho for this generation: ", best_rho)

    # select_individual will randomly choose an individual based on a roulette wheel selection
    def select_individual(self):

        individual = list

        # use random.random() to randomly choose an individual out of the population
        random_choice = self.sum_fitness * (random.random())


        # Loop through each member of the population by adding their genus values
        # once we reach a position that is >= random_choice we will stop and pick that individual
        # return the index of the individual that is chosen
        sum = 0
        i= 0
        
        for i in range(len(self.members)):
            sum = sum + self.members[i].fitness
            if(sum >= random_choice): 
                individual = i
                break
        

        return individual
    

    # mate(individual1, individual2) will perform the mating between two individuals in a population
    # The permutations will perform a "crossover" that produces two children. 
    # this method will return two permuations in the form of lists
    # this method alongside crossover() are based on Goldbergs PMX method
    def mate(self, individual_1, individual_2):

        # First we will randomly choose 2 cutoff positions in the permutations

        string_1 = self.members[individual_1].permutation
        string_2 = self.members[individual_2].permutation

        # random.sample will give 2 unique positions
        choices = random.sample(list(range(0,len(self.members[0].permutation) - 1)), 2)
        position_1 = 0
        position_2 = 0
        
        # put the two positions in order
        if(choices[0] < choices[1]):
            position_1 = choices[0]
            position_2 = choices[1]
        else:
            position_1 = choices[1]
            position_2 = choices[0]

        # We'll create two maps that save mappings from Individual1 to Individual2 and then vice versa

        # First create a substring for each individual for the area that gets singled out by the two cuts
        substring_1 = string_1[position_1:position_2 + 1]
        substring_1.pop(0)
        substring_2 = string_2[position_1:position_2 + 1]
        substring_2.pop(0)

        # map_1 will map from substr 1 to substr2
        # map_2 will map from substr2 to substr1
        map_1 = {}
        map_2 = {}

        for i in range(len(substring_1)):
            map_1[substring_1[i]] = substring_2[i]
            map_2[substring_2[i]] = substring_1[i]
        

        child_1 = self.crossover(string_1,substring_1, substring_2, map_2, position_1, position_2)
        child_2 = self.crossover(string_2, substring_2, substring_1, map_1, position_1, position_2)
        print("child 1: ", child_1)
        print("child 2: ", child_2)
    

    # the crossover function will perform the actual crossover between two permutations
    # to generate a single child from two parents
    # to generate two children from two parents, the function will be called twice within crossover()

    #takes a two lists and two dictionaries as args corresponding to first parent, 2nd parent, and their mappings
    # from each parent to the other


    def crossover(self, string, substring_1, substring_2, map, position_1, position_2):

        # extract the base of the child string from parent 
        # we'll use the map constructed from the other parent
        # to fill in the values that will get repeated after inserting the substring from that parent into
        # the middle of the base string
        child = string[:position_1 + 1] + string[position_2 + 1 :]

        seg_1 = string[:position_1 + 1]
        seg_2 = string[position_2 + 1:]

        #print("substring:", substring)

        map_keys = list(map.keys())

        # Loop through and see if the value is in list(map.keys())
        # for both segments
        for i in range(len(seg_1)):
            while(seg_1[i] in map_keys):
                seg_1[i] = map[seg_1[i]]

        for i in range(len(seg_2)):
            while(seg_2[i] in map_keys):
                seg_2[i] = map[seg_2[i]]



        # Now we will concatenate the strings together

        offspring = seg_1 + substring_2
        offspring = offspring + seg_2

        
        return offspring


    # Generate population will contain all of the steps for generating a new population
    # including selecting individuals to mate in current population, crossing over the permutations
    # and also incorporating a mutation aspect
    # Selection, crossover and mutation, will all be in separate methods within the Population class
    def generate_population(self):

        # Have a loop that will run until the size of the new population fills up to 100 individuals

        new_members = list(range(100))

        j = 0
        # while j < size of population(100)
        while j < len(new_members):

            # parent_1 and parent_2 will be the two mates
            # call select_individual to a select a random permutation
            # parent_1, and parent_2 will be index values
            # so that can be accessed from population for crossover
            parent_1 = self.select_individual()
            parent_2 = self.select_individual()
            #print(parent_1)
            #print(parent_2)

            print("Parent 1: ", self.members[parent_1].permutation)
            print("Parent 2: ", self.members[parent_2].permutation)

            self.mate(parent_1, parent_2)

            j = j + 2



        # first we need to select the next set of individuals


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
    generation_size = 2
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


    ringel_genus = optimal_genus(n)
    population.set_ringelgenus(ringel_genus)

    # calculate metrics of the initial population

    population.metrics()

    # Generate the next population of individuals
    # we'll do the genetic algorithm for the size indicated by generationsize

    for i in range(generation_size):
        # generate new population
        population.generate_population()

        # calculate metrics
        population.metrics()

        # continue on to next generation
    

    # print final results
    # population.results()
   

# MAIN

print("Input an n: ")
n = int(input())

ga_search(n)





