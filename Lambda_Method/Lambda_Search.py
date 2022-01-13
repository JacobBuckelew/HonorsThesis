"""
Lambda_Search.py takes as input a finite cyclic group of the form 12m + 7 and finds a set of lambdas
such that each lambda is composed of triplets.  

Jacob Buckelew
Rollins College
Spring 2022

"""

# lambda_search will perform a brute force style search to find a set of lambdas in the form of triplets

def lambda_search(n):

    # Save lambdas in a list

    lambdas = []

    current_lambda = []
        









# Get input from user

print("Input an n to indicate a specific finite cyclic group of the form 12n + 7:")

n = float(input())

n = int((12.0 * n) + 7)

print("*Finding an optimal set of lambdas for Z_" + str(n) + "*")

# call function lambda_search

lambda_search(n)



