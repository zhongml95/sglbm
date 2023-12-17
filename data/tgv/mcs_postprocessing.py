import numpy as np
import matplotlib.pyplot as plt

samples_number = 100
total_time_cost = 137.54+135.969+143.776+142.977
#64 881.197+887.349
#128 8835.7
nx = ny =32
u0 = 0.01

u = np.zeros((samples_number, nx, ny))
v = np.zeros((samples_number, nx, ny))
tke = np.zeros(samples_number)

dir = "/home/zhongml95/github/sglbm/data/tgv/t5/MC"+ str(nx)

tke_mean = np.zeros(samples_number)
tke_var = np.zeros(samples_number)

for n in range(0, samples_number):
    filename_u = dir+"/u_" + str(n) + ".dat"
    filename_v = dir+"/v_" + str(n) + ".dat"

    u[n, :, :] = np.loadtxt(filename_u)
    v[n, :, :] = np.loadtxt(filename_v)

    for i in range(0, nx):
        for j in range(0, ny):
            tke[n] += (u[n,i,j]*u[n,i,j] + v[n,i,j]*v[n,i,j]) * 2 / (nx*ny*u0*u0)
    
    if (n == 0):
        tke_mean[n] = tke[n]
        tke_var[n] = 0
    else:
        tke_mean[n] = np.mean(tke[0:n])
        tke_var[n] = np.var(tke[0:n])

tke_mean_err = np.abs((tke_mean - tke_mean[-1]) / tke_mean[-1])
tke_var_err = np.abs((tke_var - tke_var[-1]) / tke_var[-1])

time_list = np.linspace(total_time_cost / samples_number, total_time_cost, samples_number)

saved_filename = "tke_" + str(nx) + ".dat"
np.savetxt(saved_filename, tke)



plt.figure()
plt.plot(time_list, tke_mean_err, label=r'$nx=32, mean$')
plt.plot(time_list, tke_var_err, label=r'$nx=32, std$')

# Setting x-axis as log scale
plt.legend(fontsize = 'large')
plt.yscale('log')

plt.grid(True)

# Setting labels and title
plt.xlabel('time cost')
plt.ylabel('relative error')
plt.show()
# Display the plot
plt.savefig("./results/computational_cost.png")



