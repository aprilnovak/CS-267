from matplotlib import pyplot as plt
import numpy as np
import sys

filename = sys.argv[1]

intermediate = []
data = []
mpi = []
nel = []

serial = []
parallel = []
elems = []
procs = []
speedup = []
with open(filename, 'r') as f:
	for line in f:
		raw = line.split(' ')
		intermediate = [raw[i].rstrip() for i in range(0, len(raw)) if raw[i] != '']
		data.append(intermediate[6])
		mpi.append(intermediate[2])
		nel.append(intermediate[4])

for i in range(len(data)):
	if (i % 2 == 0):
		parallel.append(data[i])
		elems.append(nel[i])
		procs.append(mpi[i])
	else:
		serial.append(data[i])
		speedup.append(float(data[i])/float(data[i - 1]))


data = [float(data[i]) for i in range(0, len(data))]
mpi = [float(mpi[i]) for i in range(0, len(mpi))]
nel = [float(nel[i]) for i in range(0, len(nel))]
#x    = np.linspace(0.0, 1.0, len(data[0:nodes]))

#plt.plot(procs, parallel, 'o-', label='parallel')
#plt.plot(procs, serial, 'o-', label='serial')
#plt.plot(procs, speedup, label='speedup')
#plt.xlabel('Number of MPI Processes')
#plt.ylabel('Runtime (s)')
#plt.ylabel('Speedup')

# two axes
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(mpi, parallel, label='parallel')
ax1.set_ylabel('Parallel Runtime (s)')

ax2 = ax1.twinx()
ax2.plot(mpi, serial, 'r-')
ax2.set_ylabel('Serial Runtime (s)', color='r', label='serial')
for tl in ax2.get_yticklabels():
    tl.set_color('r')
plt.legend()
#plt.show()
plt.savefig('runtime-multi-node.png')
