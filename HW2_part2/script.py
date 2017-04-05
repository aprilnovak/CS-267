from matplotlib import pyplot as plt

#taskspercore = [1, 2, 3, 4, 5]
#time = [1.44465, 1.32905, 1.32644, 1.58362, 1.44324]

#plt.plot(taskspercore, time, 'o-')
#plt.xlabel('Tasks per Core for 24 Cores')

#nodes = [1, 1, 2, 4, 6, 12, 18, 24]
#time = [8.08386, 30.6676, 18.0021*2.0, 11.6832*4.0, 9.6257*6.0, 7.88603*12.0, 7.16461*18.0, 6.26946*24.0]
#time = [8.08386, 30.6676, 18.0021, 11.6832, 9.6257, 7.88603, 7.16461, 6.26946]

N = [100000, 200000, 300000, 600000, 900000, 1200000]
time = [10.172, 11.581, 12.7793, 17.7378, 22.4619, 25.5497]


plt.plot(N, time, 'o-')
plt.xlabel('Number of Particles')
plt.ylabel('Runtime (s)')
plt.savefig('weakscaling_mpidebugging.png')
