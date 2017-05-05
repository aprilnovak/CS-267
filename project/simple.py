from matplotlib import pyplot as plt
import numpy as np

# data used to show that CG scales as N^2
scalingel = [1000, 2000, 4000, 8000, 16000]
scalingt = [0.024, 0.104, 0.4125, 1.656, 12.956]
scalingt2 = []

# data used for runtime for MPI, 12 nodes
mpi = [12, 24, 36, 48, 60, 72, 84, 96]
time1 = [1.864, 0.484, 0.248, 0.148, 0.112, 0.092, 0.076, 0.068]
time2 = [500.3, 504.1, 520.15, 504, 496, 501, 530, 489]
raw = []

# data used for single-node tests
mpi = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]
time1 = [16.309, 4.228, 2.008, 1.18, 0.772, 0.54, 0.4, 0.312, 0.248, 0.208, 0.172, 0.148]
time2 = [65.52, 65.53, 69.58, 65.23, 67.34, 63.12, 61, 62, 61, 64, 65, 65.4]
raw = []
scaled = []

for i in range(len(mpi)):
	raw.append(time2[i]/time1[i])
	scaled.append(raw[i]/(mpi[i]*mpi[i]))

# weak scaling
nodes = [2, 4, 8, 16]
time1 = [10.512/10.512, 10.49/10.512, 10.52/10.512, 10.53/10.512]
time2 = [10.452/10.452, 10.532/10.452, 10.555/10.452, 10.576/10.452]	
time3 = [10.856/10.856, 10.8846/10.856, 10.8606/10.856, 10.888/10.856]

# strong scaling, 100,000 elements
nodes= [2, 4, 8, 16]
time1 = [65.53/65.53, 65.53/16.33, 65.53/4.1, 65.53/1.084]
time2 = [16.353/16.353, 16.353/4.148, 16.353/1.048, 16.353/0.3]
time3 = [7.5244/7.5244, 7.5244/1.904, 7.5244/0.5, 7.5244/0.152]
#plt.plot(nodes, time1, 'o-', label='1 task/node')
#plt.plot(nodes, time2, 'o-', label='2 task/node')
#plt.plot(nodes, time3, 'o-', label='3 task/node')
#plt.xlabel('Number of Nodes')
#plt.ylabel('Speedup')
#plt.legend()
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.plot(mpi_2, time1_2_scaled, 'o-')
#ax1.set_ylabel('Parallel Runtime (s)')
#
#ax2 = ax1.twinx()
#ax2.plot(mpi_2, time2_2, 'ro-')
#ax2.set_ylabel('Serial Runtime (s)', color='r')
#for tl in ax2.get_yticklabels():
#    tl.set_color('r')



#speedup, 20,000 elements. mpinocomm removes the communications in the DD loop to estimate percentage of time spent communicating
procs = [2, 4, 6, 8, 10, 12]
nocomm = [17.32, 8.16, 5.51, 4.3, 3.22, 2.71]
comm = [39.506, 20.049, 13.808, 10.32, 8.49, 6.43]

plt.plot(procs, comm, 'o-', label='original')
plt.plot(procs, nocomm, 'o-', label='no communication')
plt.xlabel('Number of MPI Processes')
plt.ylabel('Runtime (s)')
plt.legend()
plt.savefig('fig.png')

