"""
The aim of this programme is to take raw data from a file
containing information on all three EAs. We analyse the data 
and measure diversity in the population. This is defined as 
the population average of each individual's average euclidean 
distance to all other individuals. 

"""


import numpy as np
from matplotlib import pyplot as plt
import plotly.plotly as py
import plotly.graph_objs as go
from numba import jit


S = 30  #Number of runs
G_array = np.zeros((S))  #Array that holds the number of generations per run
N = 60  #Number of individuals
D = 11  #Dimension of individual
A = 3 #Number of algorithms


#Opening files one run at the time and noting how many generations per run there were
for s in range(0, S):

    #Opening a file and copying all results
    file_name = 'results_seed%d.txt' % s
    with open(file_name) as f:
        content = f.readlines()
    f.close()

    #Specifically checking for the number of runs
    for line in range(0, len(content)):
        if(content[line].split()[0] == '===============NewGeneration'):
            line += 1
            G_array[s] = int(content[line].split()[0])


    #Clearing some memory
    del content

#Checking whether the number of generations is the same for each run
#If not, output some text to make the user aware of this (it might lead to some discrepancies later on)
for s in range(0, S):
    if(G_array[0] == G_array[s]):
        continue
    else:
        print("Pay attention")
G = int(G_array[0])

#Declaring the 3D array that will hold all results
#First index is the fitness of individual i
#Second index is what generation
#Third index is what run
resultsCA = np.zeros((D-1, N, G, S))
resultsDE = np.zeros((D-1, N, G, S))
resultsPSO = np.zeros((D-1, N, G, S))

print("Read in data")

#Now copy all the data to the results array
for s in range(0, S):

    #Opening a file and putting the information in content
    file_name = 'results_seed%d.txt' % s
    with open(file_name) as f:
        content = f.readlines()
    f.close()

    #Go through content line by line and copy the coordinates to the appropriate indices
    line = 3
    for t in range(0, G):
        for algorithm in range(A):
            if algorithm == 0:   
                for i in range(0, N):
                    temp = content[line].strip(',').strip('[').strip(']').replace(",","").split()[:D-1]
                    resultsCA[:,i, t, s] = temp
                    line += 1
                line-=1
            if algorithm == 1:
                for i in range(0, N):
                    resultsDE[:,i, t, s] = content[line].strip(',').strip('[').strip(']').replace(",","").split()[:10]
                    line += 1
                line-=1
            if algorithm == 2:
                for i in range(0, N):
                    resultsPSO[:,i, t, s] = content[line].strip(',').strip('[').strip(']').replace(",","").split()[:10]
                    line += 1
                line+=1
            line += 2  #Skipping lines so that we skip the information on the generation number

    #Clearing some memory
    del content






#Now it's time for some statistical analysis
#First take the average per generation for each run
diversity_CA_per_gen = np.zeros((G,S))
diversity_DE_per_gen = np.zeros((G,S))
diversity_PSO_per_gen = np.zeros((G,S))

# analysis of CA

for s in range(0, S):
    for t in range(0, G):
        # mid = np.zeros(D-1)
        # for d in range(D-1):
        #     mid[d] = np.mean(resultsCA[d,:,t,s])
        distances = 0
        for ind in range(N):
            for ind2 in range(N):
                distances += np.linalg.norm(resultsCA[:,ind, t, s] - resultsCA[:,ind2, t, s])
        diversity_CA_per_gen[t, s] = distances / (N*N)

# analysis of DE
print(diversity_CA_per_gen[:,1])
for s in range(0, S):
    for t in range(0, G):
        # mid = np.zeros(D-1)
        # for d in range(D-1):
        #     mid[d] = np.mean(resultsDE[d,:,t,s])
        distances = 0
        for ind in range(N):
            for ind2 in range(N):
                distances += np.linalg.norm(resultsDE[:,ind, t, s] - resultsDE[:,ind2, t, s])
        diversity_DE_per_gen[t, s] = distances / (N*N)

# analysis of PSO

for s in range(0, S):
    for t in range(0, G):
        # mid = np.zeros(D-1)
        # for d in range(D-1):
        #     mid[d] = np.mean(resultsPSO[d,:,t,s])
        distances = 0
        for ind in range(N):
            for ind2 in range(N):
                distances += np.linalg.norm(resultsPSO[:,ind, t, s] - resultsPSO[:,ind2, t, s])
        diversity_PSO_per_gen[t, s] = distances / (N*N)

# Saving

np.savetxt('diversity_CA_per_gen.txt',diversity_CA_per_gen, fmt='%f')
np.savetxt('diversity_DE_per_gen.txt',diversity_DE_per_gen, fmt='%f')
np.savetxt('diversity_PSO_per_gen.txt',diversity_PSO_per_gen, fmt='%f')

