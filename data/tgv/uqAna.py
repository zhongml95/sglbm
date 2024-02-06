import numpy as np
import matplotlib.pyplot as plt

def reshape(oldMatrix):
    size, size = oldMatrix.shape
    step = int((size-1) / 32)
    newMatrix = np.zeros([33, 33])
    for i in range(33):
        for j in range(33):
            newMatrix[i,j] = oldMatrix[step*i,step*j]
    return newMatrix


nrlist32 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
#nqlist32 = np.array([21, 21, 21, 21, 21, 21, 21, 21, 21])
#nqlist32 = np.array([1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000])
nqlist32 = 2 * nrlist32 + 1

nrlist64 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
nqlist64 = np.array([21, 21, 21, 21, 21, 21, 21, 21, 21])
#nqlist64 = np.array([100, 100, 100, 100, 100, 100, 100, 100])

err32 = np.zeros((2, nrlist32.size-1))
err64 = np.zeros((2, nrlist64.size-1))

filenameTruth32 = "/home/zhongml95/sglbm/data/tgv/t5/ReNr"+str(nrlist32[-1])+"Nq"+str(nqlist32[-1])+"N32/final/tke.dat"
TKETruth32 = np.loadtxt(filenameTruth32)
filenameTruth64 = "/home/zhongml95/sglbm/data/tgv/t5/Nr"+str(nrlist64[-1])+"Nq"+str(nqlist64[-1])+"N64/final/tke.dat"
TKETruth64 = np.loadtxt(filenameTruth64)
#TKETruth[0] = 0.36795543345255477
for i in range(nrlist32.size-1):
    filename = "/home/zhongml95/sglbm/data/tgv/t5/ReNr"+str(nrlist64[i])+"Nq"+str(nqlist32[i])+"N32/final/tke.dat"
    TKE = np.loadtxt(filename)
    err32[0, i] = abs( TKE[0] - TKETruth32[0] ) / TKETruth32[0]
    err32[1, i] = abs( TKE[1] - TKETruth32[1] ) / TKETruth32[1]

for i in range(nrlist64.size-1):
    filename = "/home/zhongml95/sglbm/data/tgv/t5/Nr"+str(nrlist64[i])+"Nq"+str(nqlist64[i])+"N64/final/tke.dat"
    TKE = np.loadtxt(filename)
    err64[0, i] = abs( TKE[0] - TKETruth64[0] ) / TKETruth64[0]
    err64[1, i] = abs( TKE[1] - TKETruth64[1] ) / TKETruth64[1]
    
# Plotting the L2 norm
plt.figure()
plt.plot(nrlist32[:-1], err32[0,:], marker='o', label=r'$nx=32, mean$')
plt.plot(nrlist32[:-1], err32[1,:], marker='o', label=r'$nx=32, std$')
plt.plot(nrlist64[:-1], err64[0,:], marker='o', label=r'$nx=64, mean$')
plt.plot(nrlist64[:-1], err64[1,:], marker='o', label=r'$nx=64, std$')

# Setting x-axis as log scale
plt.legend(fontsize = 'large')
plt.yscale('log')

plt.grid(True)

# Setting labels and title
plt.xlabel('Polynomial Order')
plt.ylabel('relative error')
plt.show()
# Display the plot
plt.savefig("/home/zhongml95/sglbm/data/tgv/sgConvergence2.png")


#plt.figure(2)
#plt.plot(nrlist[:-1], L1Err, marker='o')

# Setting x-axis as log scale
#plt.yscale('log')

# Setting labels and title
#plt.xlabel('Polynomial Order')
#plt.ylabel('L1Err')

# Display the plot
#plt.savefig("/home/zhongml95/gPCLBM/example/data/tgv/stdL1Err.png")

