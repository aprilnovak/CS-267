from matplotlib import pyplot as plt

#plt.plot([1, 2, 4, 6, 12, 16, 18, 20], [0.98, 1.90, 3.46, 5.10, 9.59, 7.63, 8.44, 9.34])
#plt.xlabel('Threads')
#plt.ylabel('Speedup')
#plt.savefig('speedup.png')


plt.plot([1, 2, 4, 6, 12, 16, 18, 20], [0.98, 10.70, 15.11, 0.76, 0.58, 0.53, 0.47, 0.29])
plt.xlabel('Threads')
plt.ylabel('Efficiency')
plt.savefig('WeakEfficiency.png')


