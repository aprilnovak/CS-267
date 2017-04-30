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


data = [float(data[i]) for i in range(0, len(data))]
mpi = [float(mpi[i]) for i in range(0, len(mpi))]
nel = [float(nel[i]) for i in range(0, len(nel))]
#x    = np.linspace(0.0, 1.0, len(data[0:nodes]))

plt.plot(procs, parallel, label='parallel')
plt.plot(procs, serial, label='serial')
plt.show()
#plt.savefig('interface1.png')
