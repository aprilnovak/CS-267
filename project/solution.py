from matplotlib import pyplot as plt
import numpy as np

nodes = 21

raw = []
intermediate = []
data = []
with open('output.txt', 'r') as f:
	for line in f:
		raw = line.split(' ')
		intermediate = [raw[i] for i in range(0, len(raw) - 1) if raw[i] != '']
		data.extend(intermediate)

data = [float(data[i]) for i in range(0, len(data))]
x    = np.linspace(0.0, 1.0, len(data[0:nodes]))

for i in range(0, len(data), nodes):
	plt.plot(x, data[i:(i + nodes)], '*-')

plt.show()
