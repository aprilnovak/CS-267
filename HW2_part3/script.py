from matplotlib import pyplot as plt

#taskspercore = [1, 2, 3, 4, 5]
#time = [1.44465, 1.32905, 1.32644, 1.58362, 1.44324]

#plt.plot(taskspercore, time, 'o-')
#plt.xlabel('Tasks per Core for 24 Cores')

#nodes = [1, 1, 2, 4, 6, 12, 18, 24]
#time = [8.08386, 30.6676, 18.0021*2.0, 11.6832*4.0, 9.6257*6.0, 7.88603*12.0, 7.16461*18.0, 6.26946*24.0]
#time = [8.08386, 30.6676, 18.0021, 11.6832, 9.6257, 7.88603, 7.16461, 6.26946]

#N = [100000, 200000, 300000, 600000, 900000, 1200000]
#time = [10.172, 11.581, 12.7793, 17.7378, 22.4619, 25.5497]

#N = [1000, 2000, 4000, 8000, 16000, 32000]
#t1 = [0.576376, 1.10973, 2.0688, 4.01724, 8.17799, 36.0307]
#t2 = [0.78472, 0.986641, 0.993625, 0.990414, 1.2043, 1.28003]
#t3 = [0.799863, 0.959368, 0.981616, 0.986997, 1.21132, 1.28564]

N = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]
t = [0.78472, 0.986641, 0.993625, 0.990414, 1.2043, 1.28003, 1.43159, 1.98479, 3.3693, 5.62326, 10.0881, 19.0657, 37.5161, 74.3253, 147.649]

Nlog = []
for i in range(len(N)):
	Nlog.append(N[i]*0.005)

#plt.plot(N, t1, 'o-')
#plt.loglog(N, t, 'o-')
#plt.loglog(N, Nlog, 'r-')
#plt.xlabel('Number of Particles/1000')
#plt.ylabel('Runtime (s)')
#plt.savefig('gpu_final.png')

N = [1, 2, 4, 6, 9, 12, 15, 18, 21, 24]
N2 = [1, 2, 4, 6, 12, 18, 24]
sopenmp = [0.93, 1.78, 3.27, 4.59, 6.63, 8.7, 6.26, 7.24, 8.59, 10.93]
smpi = [1, 1.98, 3.89, 5.66, 10.51, 15.21, 18.34]

fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(N, sopenmp, 'g-', label='OpenMP')
ax1.plot(N2, smpi, 'b-', label='MPI')
ax2.plot([1], [816.25], 'ro', label='CUDA')

ax1.set_xlabel('Number of Threads (OpenMP) or Processes (MPI)')
ax1.set_ylabel('Speedup (OpenMP or MPI)', color='g')
ax2.set_ylabel('Speedup (CUDA)', color='r')
plt.legend()
plt.savefig('speedups.png')
