from matplotlib import pyplot as plt
import numpy as np

raw = []
with open('output.txt', 'r') as f:
	for line in f:
		raw = line.split(' ')

data = [raw[i] for i in range(0, len(raw) - 1) if raw[i] != '']
data = [float(data[i]) for i in range(0, len(data))]

x = np.linspace(0.0, 1.0, len(data))
plt.plot(x, data, '*-')
plt.show()
