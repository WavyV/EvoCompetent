import numpy as np
from matplotlib import pyplot as plt


function = 'katsuura' #'bentcigar' or 'katsuura' or 'schaffers'
algorithms = ['unimodal', 'DE', 'PSO']

averages = dict()
champions = dict()

for algorithm in algorithms:
    averages[algorithm] = np.loadtxt('averages_%s_%s.txt' % (function, algorithm), dtype=float)
    champions[algorithm] = np.loadtxt('champions_%s_%s.txt' % (function, algorithm), dtype=float)

for algorithm in algorithms:
    G = len(champions[algorithm])

    champion_median = np.zeros((G))
    champion_firstquartile = np.zeros((G))
    champion_thirdquartile = np.zeros((G))
    average_median = np.zeros((G))
    average_firstquartile = np.zeros((G))
    average_thirdquartile = np.zeros((G))
    for t in range(0, G):
        champion_median[t] = np.median(champions[algorithm][t, :])
        champion_firstquartile[t] = np.percentile(champions[algorithm][t, :], 25)
        champion_thirdquartile[t] = np.percentile(champions[algorithm][t, :], 75)
        average_median[t] = np.median(averages[algorithm][t, :])
        average_firstquartile[t] = np.percentile(averages[algorithm][t, :], 25)
        average_thirdquartile[t] = np.percentile(averages[algorithm][t, :], 75)

    if algorithm == 'unimodal':
        color = 'blue'
    if algorithm == 'DE':
        color = 'red'
    if algorithm == 'PSO':
        color = 'green'

    plt.plot(range(0,G), average_median, color=color)
    plt.fill_between(range(0, G), average_thirdquartile, average_firstquartile, color=color, alpha=0.15, linewidth=0)
    plt.plot(range(0,G), champion_median, color=color, linestyle='--')
    plt.fill_between(range(0, G), champion_firstquartile, champion_thirdquartile, color=color, alpha=0.15, linewidth=0)

plt.xlabel('Generation')
plt.ylabel('Fitness')
plt.xlim(xmin=0, xmax=10000)
plt.ylim(ymin=0, ymax=10.2)
plt.savefig('%s.png' % function, bbox_inches='tight')
plt.show()
