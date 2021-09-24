""" 

CayleyMap_Zn_HC.py applies a hill climbing approach to solving the problem of minimizing 
of a Cayley Map embedding for the finite cyclic group Zn 

Jacob Buckelew
Fall 2021
Rollins College
"""

# optimal_genus calculates the optimal genus for a complete graph Kn embedding
# Based on the Ringel and Young Theorem
# returns an int
def optimal_genus(n):
    genus = math.ceil(((n-3)*(n-4))/12)
    return genus



# hc_search will perform a local search using the hill climbing algorithm
# to find an optimal cayley map genus for the given Zn group

# will return a list representing rho and the genus of rho as an int
def hc_search(group):

    n = group_size
    


# MAIN

# Ask for user input to decide on which Zn group we are working with
print("Input an n to indicate which Z_n group you would like to calculate an optimal genus for:")
n = int(input)

# call hc_search(n) which will carry the bulk of the Hill climbing algorithm
# returns the best rho the algorithm could find and its genus

best_rho, rho_genus = hc_search(n)
