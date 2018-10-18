import numpy as np
from matplotlib import pyplot as plt


S = 100  #Number of runs
G_array = np.zeros((S))  #Array that holds the number of generations per run
N = 100  #Number of individuals
D = 11  #Dimension of individual


#Opening files one run at the time and noting how many generations per run there were
for s in range(0, S):

    #Opening a file and copying all results
    file_name = 'results_katsuura/results_seed%d.txt' % s
    with open(file_name) as f:
        content = f.readlines()
    f.close()

    #Specifically checking for the number of runs
    for line in range(0, len(content)):
        if(content[line].split()[0] == '=================='):
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
results = np.zeros((N, G, S))


#Now copy all the data to the results array
for s in range(0, S):

    #Opening a file and putting the information in content
    file_name = 'results_katsuura/results_seed%d.txt' % s
    with open(file_name) as f:
        content = f.readlines()
    f.close()

    #Go through content line by line and copy the coordinates to the appropriate indices
    line = 2
    for t in range(0, G):
        for i in range(0, N):
            results[i, t, s] = content[line].split()[0].strip(',').strip('[').strip(']')
            line += 1
        line += 2  #Skipping lines so that we skip the information on the generation number

    #Clearing some memory
    del content


#Now it's time for some statistical analysis
#First take the average per generation for each run
average_per_generation_one_run = np.zeros((G, S))
champions_per_run = np.zeros((G, S))
for s in range(0, S):
    for t in range(0, G):
        average_per_generation_one_run[t, s] = np.mean(results[:, t, s])
        champions_per_run[t, s] = np.max(results[:, t, s])

np.savetxt('results/averages_katsuura_DE.txt', average_per_generation_one_run, fmt='%f')
np.savetxt('results/champions_katsuura_DE.txt', champions_per_run, fmt='%f')


overall_champions = np.zeros((G))
for t in range(0, G):
    overall_champions[t] = np.max(champions_per_run[t, :])

median_per_generation = np.zeros((G))
first_quartile = np.zeros((G))
third_quartile = np.zeros((G))
for t in range(0, G):
    median_per_generation[t] = np.median(average_per_generation_one_run[t, :])
    first_quartile[t] = np.percentile(average_per_generation_one_run[t, :], 25)
    third_quartile[t] = np.percentile(average_per_generation_one_run[t, :], 75)

#
# np.savetxt('results/median_katsuura_DE.txt', median_per_generation, fmt='%f')
# np.savetxt('results/firstquartile_katsuura_DE.txt', first_quartile, fmt='%f')
# np.savetxt('results/thirdquartile_katsuura_DE.txt', third_quartile, fmt='%f')
# np.savetxt('results/champions_katsuura_DE.txt', overall_champions, fmt='%f')


#Plotting median with interquartile distance
plt.plot(range(0,G), median_per_generation, color='blue', label='Median')
plt.fill_between(range(0, G), third_quartile, first_quartile, color='blue', alpha=0.15, linewidth=0)
plt.plot(range(0,G), overall_champions, color='blue', linestyle='--', label='Best')
plt.xlabel('Generation')
plt.ylabel('Fitness')
plt.xlim(xmin=0)
plt.ylim(ymin=0)
plt.legend()
plt.show()
