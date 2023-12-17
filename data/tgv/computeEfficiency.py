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

nx = 32
samples_number = 10000

nrlist = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10])#, 13, 14, 15, 16, 17, 18])
nqlist = np.array([101, 101, 101, 101, 101, 101, 101, 101, 2000, 2000, 10000])#, 27, 29, 31, 33, 35, 37])
#nqlist32 = np.array([100, 100, 100, 100, 100, 100, 100, 100])

err_gpc = np.ones((2, nrlist.size))
err_mcs = np.zeros((2, samples_number))
timeCost_gpc = np.zeros(nrlist.size)
timeCost_mcs = np.zeros(samples_number)

filenameTruth = "/home/zhongml95/github/sglbm/data/tgv/t5/Nr"+str(nrlist[-1])+"Nq"+str(nqlist[-1]) +"N" + str(nx) + "/final/tke.dat"
TKETruth = np.loadtxt(filenameTruth)
#TKETruth[0] = 0.36795543345255477
for i in range(1,nrlist.size):
    filename = "/home/zhongml95/github/sglbm/data/tgv/t5/Nr"+str(nrlist[i-1])+"Nq"+str(nqlist[i-1]) +"N" + str(nx) + "/final/tke.dat"
    TKE = np.loadtxt(filename)
    err_gpc[0, i] = abs( TKE[0] - TKETruth[0] ) / TKETruth[0]
    err_gpc[1, i] = abs( TKE[1] - TKETruth[1] ) / TKETruth[1]
    timeCost_gpc[i] = TKE[3]


mcs_file_name = "/home/zhongml95/olb/apps/liang/tgv2d/tke_" + str(nx+1) + ".dat"
tke_mcs = np.loadtxt(mcs_file_name)
tke_mean_mcs = np.zeros(samples_number)
tke_var_mcs = np.zeros(samples_number)

for n in range(0, samples_number):    
    if (n == 0):
        tke_mean_mcs[n] = tke_mcs[n]
        tke_var_mcs[n] = 0
    else:
        tke_mean_mcs[n] = np.mean(tke_mcs[0:n])
        tke_var_mcs[n] = np.var(tke_mcs[0:n])

err_mcs[0,:] = np.abs((tke_mean_mcs - tke_mean_mcs[-1]) / tke_mean_mcs[-1])
err_mcs[1,:] = np.abs((tke_var_mcs - tke_var_mcs[-1]) / tke_var_mcs[-1])

if (nx == 32):
    total_time_cost = 137.54+135.969+143.776+142.977
elif (nx == 64):
    total_time_cost = (881.197+887.349)*2
elif (nx == 128):
    total_time_cost = (8835.7)*4
    
time_list = np.linspace(total_time_cost / samples_number, total_time_cost, samples_number)

# Plotting the L2 norm
plt.figure()
plt.plot(timeCost_gpc, err_gpc[0,:], label=r'$SG LBM, \bar{K}(t)$')
plt.plot(timeCost_gpc, err_gpc[1,:], label=r'$SG LBM, \sigma(k(t))$')

plt.plot(time_list[0:samples_number-1], err_mcs[0,0:samples_number-1], label=r'$MCS(OpenLB), \bar{K}(t)$')
plt.plot(time_list[0:samples_number-1], err_mcs[1,0:samples_number-1], label=r'$MCS(OpenLB), \sigma(k(t))$')


# Setting x-axis as log scale
plt.legend( loc='lower center', ncol = 2, bbox_to_anchor=(0.65, 0))
#plt.xscale('log')
plt.yscale('log')

plt.grid(True)

# Setting labels and title
plt.xlabel('time cost')
plt.ylabel('relative error')
#plt.show()
# Display the plot
plt.savefig("/home/zhongml95/github/sglbm/data/tgv/computational_cost_" + str(nx) + ".png")


#plt.figure(2)
#plt.plot(nrlist[:-1], L1Err, marker='o')

# Setting x-axis as log scale
#plt.yscale('log')

# Setting labels and title
#plt.xlabel('Polynomial Order')
#plt.ylabel('L1Err')

# Display the plot
#plt.savefig("/home/zhongml95/gPCLBM/example/data/tgv/stdL1Err.png")

