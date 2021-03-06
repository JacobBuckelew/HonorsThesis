"""
Lambda_GA.py finds an optimal cayley map embedding for a given Z_n group using a genetic algorithm.


Jacob Buckelew
Spring 2022
Rollins College

"""


import math, random, copy, os, timeit, matplotlib.pyplot as plt
from os import path



# Object Classes


# Define a population class to encapsulate all of the generations in a population
# as well as best_genus for each generation, average_genus for each generation, as
# well as the most optimal running genus and its permutation throughout all generations
class Population:

    best_fitness = []
    opt_fitness = 100000
    opt_lambda = []
    avg_fitness = []
    members = list
    ringelgenus = 0
    counter = 0
    current_lambda = []
    sum_scale_fitness = 0
    fitness_mappings = {}
    group_size = 0


    def __init__(self, members, current_lambda, group_size):
        self.members = members
        self.current_lambda = current_lambda
        self.group_size = group_size
    
    def set_members(self, members):
        self.members = members

    def set_current_lambda(self, current_lambda):
        self.current_lambda = current_lambda

    def set_ringelgenus(self, ringel_genus):
        self.ringelgenus = ringel_genus

    def set_sum_scale_fitness(self, sum_scale_fitness):
        self.sum_scale_fitness = sum_scale_fitness

    def set_opt_lambda(self, opt_lambda):
        self.opt_lambda = opt_lambda

    def set_best_fitness(self, best_fitness):
        self.best_fitness = best_fitness
    
    def set_opt_fitness(self, opt_fitness):
        self.opt_fitness = opt_fitness

    def set_avg_fitness(self, avg_fitness):
        self.avg_fitness = avg_fitness

    def set_fitness_mappings(self, fitness_mappings):
        self.fitness_mappings = fitness_mappings
    

    def update_counter(self, n):
        self.counter = self.counter + n
    

    
    # metrics will iterate through the population to find avg_genus, best_genus, opt_genus, and opt_rho
    # for each generation
    def metrics(self, solution_start):


        best_fitness = 100000

        total_fitness = 0
        # iterate through the members
        count = 0
        for member in self.members:
            # find individual's genus
            fitness = member.calculate_fitness(member.lam_orient, member.lam_map, member.lam)
            member.set_fitness(fitness)
            

            
            # If statement to determine best genus so far
            # set best_genus if member fitness is less than running best genus for this generation
            if(member.fitness < best_fitness):
                # update population's optimal rho and genus
                best_lambda = member.lam
                best_fitness = member.fitness
            
            # Print if we see a rho that satisfies ringel definition, this will be most optimal genus and rho
            # or if it is less than overall optimal genus found by the algorithm
            if(member.fitness <= self.opt_fitness):
                if(member.fitness == 1 and (member.lam not in self.opt_lambda)):
                    #print("Found Ringel Optimal Solution")
                    #print("Lambda:", member.lam)
                    #print("\n")
                    #print("Cycles in Rho:", member.fitness)
                    #print("\n")
                    # lam = member.lam
                    # print(lam)
                    # print("lam:", type(lam))
                    # for i in range(len(lam)):
                    #     for j in range(len(lam[i])):
                    #         if(lam[i][j] > ((self.group_size - 1)/2)):
                    #             lam[i][j] = lam[i][j] - self.group_size
                                
                    # print("updated lam:", lam)
                    # member.set_lam(lam)
                    self.opt_lambda.append(member.lam)
                    if(len(self.opt_lambda) == 1):
                        solution_end = timeit.default_timer()
                        solution_time = solution_end - solution_start
                        print("First solution found in " + str(solution_time) + " seconds")
                        print("Lambda:", member.lam)
                    count = count + 1
                #else:
                    #print("Found more optimal lambda:")
                    #print("Lambda: ", member.lam)
                    #print("Cycles: ", member.fitness)
                    #print("\n")
                self.set_opt_fitness(fitness)

            total_fitness = total_fitness + fitness



            # add up genus to running total to find the average genus in the generation
        # print("number of optimal individuals in generation:", count)
        avg_fit = math.floor(int(total_fitness/len(self.members)))
        self.best_fitness.append(best_fitness)
        self.avg_fitness.append(avg_fit)
    
    def analyze(self, n, k):

        
        x_axis = []
        for i in range(len(self.best_fitness)):
            x_axis.append(i)
        plt.figure()
        plt.plot(x_axis, self.avg_fitness, label = "Average Fitness")
        plt.plot(x_axis, self.best_fitness, label = "Best Fitness")
        plt.xlabel("Generation")
        plt.ylabel("Cycles in Rho")
        plt.legend(bbox_to_anchor = (0.75, 1.15), ncol = 2)
        plot_dir = "Lambda_Plots/Z" + str(k)
        if(not path.exists(plot_dir)):
            print("Making New Directory")
            os.mkdir(plot_dir)    
        plt.savefig("Lambda_Plots/Z" + str(k) +"/Z" + str(k) + "Results_Lambda" + str(n) + ".jpeg")
            






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
        if(self.sum_scale_fitness == 0):
            candidates = {}

            sum_fitnesses = 0

            for i in range(len(self.members)):
                fitness = self.members[i].fitness
                scale_fitness = math.floor((1.0/fitness) * (self.group_size))
                if(scale_fitness not in list(candidates.keys())):
                    candidates[scale_fitness] = []
                    sum_fitnesses = sum_fitnesses + scale_fitness
                
                candidates[scale_fitness].append(self.members[i])
            
            self.set_fitness_mappings(candidates)
            self.set_sum_scale_fitness(sum_fitnesses)

            
            
        # # use random.random() to randomly choose a genus out of the population
        random_choice = self.sum_scale_fitness * (random.random())

        # # search through keys

        keys = list(self.fitness_mappings.keys())
        sum = 0
        candidate = Individual

        for i in range(len(keys)):
            sum = sum + keys[i]
            if(sum >= random_choice):
                # choose the genus corresponding to the value that was added to the sum
                pool = self.fitness_mappings[keys[i]]
                random_individual = random.randint(0, len(pool) - 1)
                candidate = pool[random_individual]
                break



        return candidate
    

    def save_individuals(self):
    # We'll save the top 10% of individuals
    # Use python's sorted function to sort objects by their fitness values
        k = math.floor(.10 * len(self.members))
        
        individuals = sorted(self.members, key =lambda individual: individual.fitness)
        top_k = list
        if(individuals[k].fitness == self.best_fitness[self.counter - 1]):
            top_individuals = []
            for i in range(len(individuals)):
                if(individuals[i].fitness != self.best_fitness[self.counter - 1]):
                    break
                top_individuals.append(individuals[i])
            top_k = random.sample(top_individuals, k)
        else:
            top_k = individuals[0:k]

        return top_k


    # mate(individual1, individual2) will perform the mating between two individuals in a population
    # The permutations will perform a "crossover" that produces two children. 
    # this method will return two permuations in the form of lists
    # this method alongside crossover() are based on Goldbergs PMX method
    def mate(self, individual_1, individual_2):

        # First we will randomly choose 1 cutoff positions in the permutations

        string_1 = individual_1.lam_orient
        string_2 = individual_2.lam_orient

        #cutoff = math.floor(len(string_1)/2.0)

        cutoff = random.randint(0, len(string_1) - 1)

        string_1_left = string_1[0:cutoff]
        string_1_right = string_1[cutoff:]

        string_2_left = string_2[0:cutoff]
        string_2_right = string_2[cutoff:]

        child_1 = string_1_left + string_2_right
        child_2 = string_2_left + string_1_right


        return child_1, child_2
    



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

        for i in range(len(top_ten)):
            lam, lam_map = set_map(self.group_size, top_ten[i].lam_orient, self.current_lambda)
            top_ten[i].set_map(lam_map)
            top_ten[i].set_lam(lam)

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
            while(parent_2 == parent_1):
                parent_2 = self.select_individual()

            # Mate the two parents

            offspring_1, offspring_2 = self.mate(parent_1, parent_2)

            child_1 = Individual(offspring_1)

            child_2 = Individual(offspring_2)

            child_1.mutate()
            child_2.mutate()

            #mutation_p = random.random()
            #if(mutation_p < 0.01):
                #mutations = mutations + 1
                #child_1.mutate()
            
            #mutation_p = random.random()
            #if(mutation_p < 0.01):
                #mutations = mutations + 1
                #child_2.mutate()

            lam, lam_map = set_map(self.group_size, child_1.lam_orient, self.current_lambda)

            child_1.set_map(lam_map)
            child_1.set_lam(lam)

            lam, lam_map = set_map(self.group_size, child_2.lam_orient, self.current_lambda)

            child_2.set_map(lam_map)
            child_2.set_lam(lam)
    
            new_members.append(child_1)
            new_members.append(child_2)




            j = j + 2
            #print("j", j)

        # reset the members list in population
        self.members = new_members
        self.set_fitness_mappings({})
        self.set_sum_scale_fitness(0)
        #print(str(mutations) + " total mutations")

    # results() will print out the best running permutation
    # 


# Define an individual class
class Individual:
    def __init__(self, lam_orient):
        self.lam_orient = lam_orient
    
    def set_fitness(self, cycles):
        self.fitness = cycles

    def set_map(self, lam_map):
        self.lam_map = lam_map
    
    def set_lam(self, lam):
        self.lam = lam
    


    def get_fitness(self):
        return self.fitness

    
    
    
    def mutate(self):
        # Check bit by bit, and flip a bit with a small probability

        # use random.sample() to get random index
        #index = random.randint(0, len(self.lam_orient) - 1)

        p = .001

        for i in range(len(self.lam_orient)):
            random_value = random.random()
            if(random_value <= p):
                self.lam_orient[i] = (self.lam_orient[i] + 1) % 2

        
    def calculate_fitness(self, lam_orient, lam_map, lam):

        fitness = 0

        rho = []


        n = (12 * len(list(lam_map.keys()))) + 7

        
        # find rho

        for i in range(1, n):

            if(i == 1):
                cycle = [1]
            if(not self.in_rho(rho, i) and (not self.in_rho(cycle, i))):
                cycle = [i]

                for j in range(0, n):

                    next_value = lam_map[cycle[j]]
                    lam_map[cycle[j]] = 0

                    if(self.in_rho(rho, next_value) or self.in_rho(cycle, next_value)):
                        break
                    cycle.append(next_value)
            
                rho.append(cycle)
                cycle = []
            
            if(sum(lam_map.values()) == 0):
                break
        
        fitness = len(rho)
        return fitness
        
                

    def in_rho(self, rho, m):
        if(len(rho) == 0):
            return False
        if(type(rho[0]) is int):
            if(m in rho):
                return True
            else:
                return False
        else:
            for i in range(0,len(rho)):
                if(m in rho[i]):
                    return True
            return False



# Functions

# optimal_genus calculates the optimal genus for a complete graph Kn embedding
# Based on the Ringel and Young Theorem
# returns an int
def optimal_genus(n):
    genus = math.ceil(((n-3)*(n-4))/12)
    return genus


def random_orient(n):

    triples = int((n-1)/3)

    random_num = random.random()
    orientation = []
    for i in range(0,triples):
        if(random_num < 0.5):
            orientation.append(0)
        else:
            orientation.append(1)
        random_num = random.random()

    
    return orientation




def in_lambda(lam, m):
    if len(lam) == 0:
        return False
    
    return (m in lam[0]) or in_lambda(lam[1:], m)


def set_map(n, orient, current_lambda):


    triples = int((n-1)/3)

    lam_map = {}

    lam = copy.deepcopy(current_lambda)


    for i in range(0, len(lam)):
        if(orient[i] == 0):
            continue
        else:
            lam[i].reverse()

    if(len(lam) == triples):
        for i in range(0, triples):
            for j in range(0, 3):
                i_inv = n -  lam[i][j]
                if i_inv in list(lam_map.keys()):
                    continue
                # edge case
                if(j == 2):
                    lam_map[i_inv] = lam[i][0]
                else:
                    lam_map[i_inv] = lam[i][j + 1]

    return lam, lam_map



# set_lambda will perform a brute force style search to find a set of lambdas in the form of triplets

def set_lambda(n, lam):

    triples = int((n - 1)/3)




    # initialize a new lambda if starting the program fresh
    if(len(lam) == 0):
        current_lambda = [[1, 2, n - 3]]
        x = 1
        y = 2
    else:
        current_lambda = lam
        
        # Backtrack to previous factor in lambda before beginning search

        current_lambda.pop()
        temp = current_lambda.pop()

        x = temp[0]
        y = temp[1] + 1

    if(len(current_lambda) == 0):
        return current_lambda

    start = timeit.default_timer()
    while(current_lambda[0][1] < current_lambda[0][2]):
        # first check if lambda is already full, in which case we have a full set of triplets
        if len(current_lambda) == triples:
            break
        
        # Now for the search process

        for i in range(x, n):
            # if the i is not in lambda, we'll consider the i for this factor in lambda
            if not in_lambda(current_lambda, i):
                # now find which j we will use 

                for j in range(y,n):
                    if not in_lambda(current_lambda, j) and not in_lambda(current_lambda, i) and j > i:

                        # figure out the missing k required to get a triplet

                        k = (n - j - i) % n


                        # if k has not been used we have a valid triplet and then we append the factor to lambda
                        if (k > j) and (not in_lambda(current_lambda, k)) and (not in_lambda(current_lambda, i)) and (not in_lambda(current_lambda, j)):

                            current_lambda.append([i, j, k])   
            # Need to check if the i was actually added to lambda 
            # There is a possibility that the factor did not append 

            if not in_lambda(current_lambda, i):

                temp = current_lambda.pop()

                x = temp[0]

                y = temp[1] + 1


                # if resetting the factors in lambda, we need to create the first factor now

                if len(current_lambda) == 0:

                    z = (n - x - y) % n

                    current_lambda.append([x,y,z])

                    x = x + 1

                    y = x + 1
                break

    end = timeit.default_timer()
    print("Time elapsed to find lambda: " + str(end - start) + " seconds\n")

    return current_lambda

def results(n, results_map):
    
    print("***********************************************\n")
    print("RESULTS FOR Z_" + str(n) +"\n***********************************************\n")

    keys = results_map.keys()


    for key in keys:

    
        print("Lambda: " + str(key))
        print("\n")
        # print("Optimal Orientations(s):\n")
        if(len(results_map[key][0]) == 0):
            print("No optimal orientations.\n")
            print("Least number of cycles generated: " + str(results_map[key][1]))
            print("\n")
            print("----------------------------------------------\n")
        else:
            #print("Lambdas: ", results_map[key][0])
            print("Number of optimal orientations: " + str(len(results_map[key][0])) + "\n")
            print("One optimal orientation:")
            print(results_map[key][0])
            print("\n")
            print("-----------------------------------------------")
        

def export_file(n, results_map):

    keys = results_map.keys()
    
    filename = "Lambda_Results/LambdaResults_Z" + str(n) + ".txt"

    #if os.path.exists(filename):
        #append_write = 'a' # append if already exists
    #else:
        #append_write = 'w' # make a new file if not

    fo = open(filename, 'w')

    fo.write("***********************************************\n")
    fo.write("RESULTS FOR Z_" + str(n) +"\n***********************************************\n")

    for key in keys:


        fo.write("Lambda: " + str(key))
        fo.write("\n")
        #fo.write("Optimal Orientations(s):\n")
        if(len(results_map[key][0]) == 0):
            fo.write("No optimal orientations.\n")
            fo.write("Least number of cycles generated: " + str(results_map[key][1]))
            fo.write("\n")
            fo.write("----------------------------------------------\n")
        else:
            #fo.write("Lambdas: ", results_map[key][0])
            fo.write("Number of optimal orientations: " + str(len(results_map[key][0])) + "\n")
            fo.write("One optimal orientation:")
            fo.write(str(results_map[key][0]))
            fo.write("\n")
            fo.write("-----------------------------------------------\n")
        
    fo.write("\n\n\n")
    fo.close()



def ga_search(n, generations, individuals):

    population_size = individuals
    generation_size = generations - 1

    
    print("Would you like to save convergence graphs for each generation?")
    print("Type 1 if you would like to save graphs or else type 0.")
    response = int(input())

    solution_start = timeit.default_timer()
    # population_members will store individuals

    results_map = {}

    # keep track of average, best, most optimal running rho and its genus
    # population_members will store individuals
    population_members = []

    # for each generation set a lambda to be analyzed
    current_lambda = []

    current_lambda = set_lambda(n, current_lambda)
    #print(type(current_lambda))
    opt_fit = 100000

    simulation = 0
    while(len(current_lambda) > 1):
        simulation = simulation + 1
        #print("Initializing:")
        start = timeit.default_timer()

        for i in range(population_size):
            # initialize individuals
            # Generate random sequence of 0s and 1s which is length of lambda's size

            orientation = random_orient(n)

            member = Individual(orientation)

            lam, lam_map = set_map(n, orientation, current_lambda)
            
            # Use the generate_lambda to initialize the population for the first generation
            member.set_map(lam_map)
            member.set_lam(lam)


            # add individual to the population of individuals
            population_members.append(member)
        
        population = Population(population_members, current_lambda, n)
        population.set_opt_fitness(opt_fit)
        population.set_opt_lambda([])
        population.set_best_fitness([])
        population.set_avg_fitness([])
        
        population.set_members(population_members)


        ringel_genus = optimal_genus(n)
        population.set_ringelgenus(ringel_genus)

        # calculate metrics of the initial population

        population.metrics(solution_start)

        population.update_counter(1)

        # Generate the next population of individuals
        # we'll do the genetic algorithm for the size indicated by generationsize

        for i in range(generation_size):
            # generate new population
            population.generate_population(population_size)

            # calculate metrics
            population.metrics(solution_start)

            # loop back and continue on to next generation
            population.update_counter(1)

        print("Num of optimal orientations:", len(population.opt_lambda))
        if(response == 1):
            population.analyze(simulation, n)
        
        opt_lam = copy.deepcopy(population.opt_lambda)
        if(len(opt_lam) > 0):
            opt_lam = opt_lam[1]
        altered_lambda = copy.deepcopy(current_lambda)
        for i in range(len(altered_lambda)):
            for j in range(len(altered_lambda[i])):
                if(altered_lambda[i][j] > ((n-1)/2)):
                    altered_lambda[i][j] = altered_lambda[i][j] - n
                if(len(opt_lam) > 0):
                    if(opt_lam[i][j] > ((n-1)/2)):
                        opt_lam[i][j] = opt_lam[i][j] - n 

        lam_tup = tuple(tuple(factor) for factor in altered_lambda)
        results_map[lam_tup] = [opt_lam, population.opt_fitness]
        current_lambda = set_lambda(n, current_lambda)
        population_members = []
    
    end_time = timeit.default_timer()
    final_time = end_time - solution_start
    print("Algorithm ran for " + str(final_time) + " seconds.")
    print("Would you like results printed out or exported?")
    print("Type 1 for printed results or 2 for exported results")
    response = int(input())
    if(response == 1):
        results(n, results_map)
    if(response == 2):
        export_file(n, results_map)

# MAIN

print("Input an n: ")
n = int(input())

n = (12 * n) + 7

print("Input a number of generations: ")

gen = int(input())

print("Input a number of individuals: ")

individ = int(input())


ga_search(n, gen, individ)
