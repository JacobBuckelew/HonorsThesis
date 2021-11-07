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

    opt_genus = 200
    best_genus = []
    avg_genus = []
    best_rho = []
    sum_scale_genus = 0
    members = list
    ringelgenus = 0
    counter = 0
    genus_mapping = {}


    def __init__(self, members):
        self.members = members
    
    def set_members(self, members):
        self.members = members

    def set_ringelgenus(self, ringel_genus):
        self.ringelgenus = ringel_genus

    def set_best_genus(self, best_genus):
        self.best_genus.append(best_genus)
    
    def set_opt_genus(self, opt_genus):
        self.opt_genus = opt_genus
    
    def set_best_rho(self, best_rho):
        self.best_genus.append(best_rho)
    
    def set_avg_genus(self, avg_genus):
        self.avg_genus.append(avg_genus)

    def set_genus_mapping(self, map):
        self.genus_mapping = map

    def set_sum_scale_genus(self, sum_genus):
        self.sum_scale_genus = sum_genus

    def update_counter(self, n):
        self.counter = self.counter + n
    
    
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
            # or if it is less than overall optimal genus found by the algorithm
            if((member.fitness == self.ringelgenus) or (member.fitness < self.opt_genus)):
                #print("Ringel and Young optimal genus found")
                #print("Genus: ", member.fitness)
                #print("Rho: ", member.permutation)
                best_genus = genus
                best_rho = member.permutation
                self.set_opt_genus(genus)

            

            # add up genus to running total to find the average genus in the generation

            avg_genus = avg_genus + genus
            
        # self.sum_scale_genus = avg_genus

        avg_genus = int(round(avg_genus/len(self.members)))
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

        # We'll scale the genus values by subtracting them from 2 times the avg genus from the previous generation
        # and then squaring the genus

        # This technique will favor the chromosomes with smaller genus values
        # We'll model selection based on a roulette wheel that chooses a particular *scaled* genus value
        # From there we can randomly choose an individual from the pool of individuals that share that genus value

        # First we'll scale and then add individuals to a map that saves them to their corresponding genus values

        # if we are running this method for the first time in this generation then we will create mapping
        # and set the sum of the scaled genuses to prevent having to do it again for selecting the other individuals
        
        if(self.sum_scale_genus == 0):
            candidates = {}

            sum_genuses = 0

            for i in range(len(self.members)):
                genus = self.members[i].fitness
                #print(self.members[i].permutation)
                #print("genus: ", genus)
                scale_genus = (((2* self.avg_genus[self.counter - 1]) - genus) ** 2)

                #print("scaled genus: ", scale_genus)

                if(scale_genus not in list(candidates.keys())):
                    candidates[scale_genus] = []
                    sum_genuses = sum_genuses + scale_genus
                
                candidates[scale_genus].append(self.members[i])
            
            self.set_genus_mapping(candidates)
            self.set_sum_scale_genus(sum_genuses)
            #print("full map: ", self.genus_mapping)

            
            
        # # use random.random() to randomly choose a genus out of the population
        random_choice = self.sum_scale_genus * (random.random())

        # # search through keys

        keys = list(self.genus_mapping.keys())
        sum = 0
        candidate = Individual

        for i in range(len(keys)):
            sum = sum + keys[i]
            if(sum >= random_choice):
                # choose the genus corresponding to the value that was added to the sum
                pool = self.genus_mapping[keys[i]]
                #print("pool:", pool)
                random_individual = random.randint(0, len(pool) - 1)
                candidate = pool[random_individual]
                break

        # print(candidate)
        # print(candidate.permutation)
        # print(candidate.fitness)

        return candidate




    def save_individuals(self):

        # We'll save the top 10% of individuals

        # Use python's sorted function to sort objects by their fitness values

        individuals = sorted(self.members, key =lambda individual: individual.fitness)
        top_ten = list
        if(individuals[10].fitness == self.best_genus[self.counter - 1]):
            top_individuals = []
            for i in range(len(individuals)):
                if(individuals[i].fitness != self.best_genus[self.counter - 1]):
                    break
                top_individuals.append(individuals[i])
                # print(i)
            top_ten = random.sample(top_individuals, 10)
        else:
            top_ten = individuals[0:10]


        #for individual in top_ten:
            #print(individual.fitness)

        return top_ten


    # mate(individual1, individual2) will perform the mating between two individuals in a population
    # The permutations will perform a "crossover" that produces two children. 
    # this method will return two permuations in the form of lists
    # this method alongside crossover() are based on Goldbergs PMX method
    def mate(self, individual_1, individual_2):

        # First we will randomly choose 2 cutoff positions in the permutations

        string_1 = individual_1.permutation
        string_2 = individual_2.permutation

        # random.sample will give 2 unique positions
        choices = random.sample(list(range(0,len(self.members[0].permutation))), 2)
        position_1 = 0
        position_2 = 0

        
        # put the two positions in order
        if(choices[0] < choices[1]):
            position_1 = choices[0]
            position_2 = choices[1]
        else:
            position_1 = choices[1]
            position_2 = choices[0]

        #print("position_1", position_1)
        #print("position_2", position_2)
  

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

        return child_1, child_2
    

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
    def generate_population(self, size):

        # Have a loop that will run until the size of the new population fills up to 100 individuals

        new_members = []

        # we'll save the top 10% of individuals from the current generation to stay in the next generation
        # in order to preserve "good" traits

        top_ten = self.save_individuals()

        new_members = new_members + top_ten


        reserved = len(new_members)

        j = 0
        mutations = 0
        # while j < size of population(100)
        while j < (size - reserved):
            # parent_1 and parent_2 will be the two mates
            # call select_individual to a select a random permutation
            # 
            parent_1 = self.select_individual()

            parent_2 = self.select_individual()
            while(parent_2.permutation == parent_1.permutation ):
                #print("Equal")
                parent_2 = self.select_individual()

            #print("Parent 1: ", self.members[parent_1].permutation)
            #print("Parent 2: ", self.members[parent_2].permutation)

            # Mate the two parents

            offspring_1, offspring_2 = self.mate(parent_1, parent_2)

            child_1 = Individual(offspring_1)

            child_2 = Individual(offspring_2)

            #print("child 1: ", child_1.permutation)
            #print("child 2:", child_2.permutation)

            # Decide whether or not they will mutate
            # if random # < 1/ length of perm size

            mutation_p = random.random()

            # set a mutation probability 
            # Use .01 to start

            if(mutation_p < 0.01):
                mutations = mutations + 1
                child_1.mutate()

            mutation_p = random.random()
            if(mutation_p < 0.01):
                mutations = mutations + 1
                child_2.mutate()
    
            new_members.append(child_1)
            new_members.append(child_2)

            print("parent1")
            print(parent_1.permutation)
            print("parent2")
            print(parent_2.permutation)
            print("child1")
            print(child_1.permutation)
            print("child2")
            print(child_2.permutation)

            child_1.set_fitness(child_1.calculate_genus(child_1.permutation))
            if(child_1.fitness < self.best_genus[self.counter - 1]):
                print("child 1 has a best genus")
                print(child_1.permutation)
                print(child_1.fitness)
                print(parent_1.permutation)
                print(parent_1.fitness)
                print(parent_2.permutation)
                print(parent_2.fitness)

            
            child_2.set_fitness(child_2.calculate_genus(child_2.permutation))


            if(child_2.fitness < self.best_genus[self.counter - 1]):
                print("child 2 has a best genus")
                print(child_2.permutation)
                print(child_2.fitness)
                print(parent_2.permutation)
                print(parent_2.fitness)
                print(parent_1.permutation)
                print(parent_1.fitness)
            




            #print(child_1.permutation)
            #print(child_2.permutation)

            if(child_1.permutation == [1,6,15,7,17,16,5,18,2,9,13,14,11,4,10,12,8,3]):
                print("Found optimal")
            if(child_2.permutation == [1,6,15,7,17,16,5,18,2,9,13,14,11,4,10,12,8,3]):
                print("Found optimal")

            j = j + 2
            #print("j", j)

        # reset the members list in population
        self.members = new_members
        self.set_genus_mapping({})
        self.set_sum_scale_genus(0)
        #print(str(mutations) + " total mutations")

    # results() will print out the best running permutation
    # 
    def results(self):

        print("Genetic Algorithm has finished")
        print("------------------------------------")
        print("Optimal rho found was " + str(self.best_rho) + " with genus " + str(self.opt_genus))
        print("Average genus values: ", self.avg_genus)
        print("Best Genus values:", self.best_genus)
        print("Final Generation:")

        for i in range(len(self.members)):
            print("Individual", self.members[i])
            print("Rho", self.members[i].permutation)
            print("Genus", self.members[i].fitness)
            print("******************************************")


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
    
    def mutate(self):
        # Choose two random indices to swap

        # use random.sample() to get unique values from 
        # 1 to len - 1

        indices = random.sample(list(range(1,len(self.permutation))), 2)

        temp = self.permutation[indices[0]]
        self.permutation[indices[0]] = self.permutation[indices[1]]
        self.permutation[indices[1]] = temp


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
        print("Found optimal")
    return permutation

    



def ga_search(n):

    population_size = 100
    generation_size = 10000
    # population_members will store individuals

    # keep track of average, best, most optimal running rho and its genus
    # population_members will store individuals
    population_members = []


    for i in range(population_size):
        # initialize individuals
        member = Individual(random_permutation(n))

        #print(member.permutation)

        # add individual to the population of individuals
        population_members.append(member)
    
    population = Population(population_members)
    
    population.set_members(population_members)


    ringel_genus = optimal_genus(n)
    population.set_ringelgenus(ringel_genus)

    # calculate metrics of the initial population

    population.metrics()

    population.update_counter(1)

    # Generate the next population of individuals
    # we'll do the genetic algorithm for the size indicated by generationsize

    for i in range(generation_size):
        # generate new population
        population.generate_population(population_size)

        # calculate metrics
        population.metrics()

        # loop back and continue on to next generation
        population.update_counter(1)
    

    # print final results
    population.results()
   

# MAIN

print("Input an n: ")
n = int(input())

ga_search(n)





