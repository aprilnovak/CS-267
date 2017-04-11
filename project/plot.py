import sys # For: command line arguments
import matplotlib # For: plotting
from matplotlib import pyplot as plt

# get the filenames and descriptions for comparison
num_files = (len(sys.argv) - 1)/2

filenames = []
descriptions = []
for i in range(0, num_files + 1, 2):
	filenames.append(sys.argv[i + 1])
	descriptions.append(sys.argv[i + 2])

my_dict = dict()
for index, filename in enumerate(filenames):
	with open(filename, 'r') as f:
		size = []
		CG = []
		total = []
		CGiterations = []
		description = descriptions[index]
		my_dict[description] = {}
		for line in f:
			data = line.split()
			print(data)
			size.append(float(data[0]))
			total.append(float(data[1]))
			CG.append(float(data[2]))
			CGiterations.append(float(data[2]) / float(data[3]))
		my_dict[description] = {'size':size, 'CG':CG, 'total':total, 'CGiterations':CGiterations}

for description in descriptions:
    plt.plot(my_dict[description]['size'], my_dict[description]['CGiterations'], label=description)
    #plt.plot(my_dict[description]['size'], my_dict[description]['CG'], label='CG solver')
    plt.ylabel('Runtime per CG iteration (s)')
    plt.xlabel('Size')
plt.legend()
plt.savefig(description + '.png')
