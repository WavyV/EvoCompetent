import numpy as np
import matplotlib.pyplot as plt

G = 665 #Number of Generations
S = 5 #number of seeds used in plotting
A = 3 #Number of algorithms

algorithms = ['CA', 'DE', 'PSO']

# Read and store data

diversity_per_gen = np.zeros((len(algorithms),G, S))

for p, a in enumerate(algorithms):
	file_name = 'diversity_%s_per_gen.txt' %a
	with open(file_name) as f:
		content = f.readlines()
	f.close()

	for i, line in enumerate(content):
		for j, token in enumerate(line.split(" ")):
			diversity_per_gen[p, i, j] = float(token)

# compute means and standard deviations
e = np.zeros((A,G))
mean = np.zeros((A, G))
for a in range(len(algorithms)):
	for gen in range(G):
		e[a, gen] = np.std(diversity_per_gen[a,gen,:])
		mean[a, gen] = np.mean(diversity_per_gen[a,gen,:])


for p, a in enumerate(algorithms):
	x = range(G)
	#plt.errorbar(x, mean[p,:], e[p,:])
	plt.plot(x, mean[p,:], label="%s"%a)

plt.yscale('linear')	
plt.title("Diversity per generation in %s" %a)
plt.xlabel("Generations")
plt.ylabel("Diversity")
plt.legend()
plt.show()

