"""
Lambda_Search.py takes as input a finite cyclic group of the form 12m + 7 and finds a set of lambdas
such that each lambda is composed of triplets. Code in inspired by Miriam Scheinblum's brute force algorithm
used in her Honors Thesis titled Optimal Cayley Map Embeddings of Complete Graphs.

Jacob Buckelew
Rollins College
Spring 2022

"""

def in_lambda(lam, m):
    if len(lam) == 0:
        return False
    
    return (m in lam[0]) or in_lambda(lam[1:], m)



# lambda_search will perform a brute force style search to find a set of lambdas in the form of triplets

def lambda_search(n, lam):

    # Save lambdas in a list

    triples = int((n - 1)/3)


    # initialize a new lambda if starting the program fresh
    if(len(lam) == 0):
        current_lambda = [[1,2, n - 3]]
        x = 1
        y = 2
    else:
        current_lambda = lam
        
        # Backtrack to previous factor in lambda before beginning search

        current_lambda.pop()
        temp = current_lambda.pop()

        x = temp[0]
        y = temp[1] + 1



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
    
    return current_lambda
    


# Get input from user

print("Input an n to indicate a specific finite cyclic group of the form 12n + 7:")

n = float(input())

n = int((12.0 * n) + 7)

print("*Finding an optimal set of lambdas for Z_" + str(n) + "*")

# call function lambda_search

lam = []

lam = lambda_search(n, lam)

while(len(lam) != 1):
    print(lam)

    lam = lambda_search(n, lam)






