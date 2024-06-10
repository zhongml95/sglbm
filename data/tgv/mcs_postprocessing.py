import numpy as np
import matplotlib.pyplot as plt
import os

def save_to_tecplot(resolution,u_mean, v_mean, u_std, v_std, filename, title="Velocity Field"):
    """
    Save data to a Tecplot format .dat file.

    Parameters:
    - resolution: resolution of the grid in the x or y direction.
    - u_mean: 2D array of mean velocity in horizontal direction.
    - v_mean: 2D array of mean velocity in vertical direction.
    - u_std: 2D array of std velocity in horizontal direction.
    - u_std: 2D array of std velocity in vertical direction.
    - filename: Path and name of the file to save.
    - title: Title of the dataset (optional).
    - variable_name: Name of the variable for the QoI (optional).
    """
    with open(filename, 'w') as f:
        f.write(f"TITLE = \"{title}\"\n")
        f.write(f"VARIABLES = \"X\" \"Y\"  \"umean\" \"vmean\" \"ustd\" \"vstd\" \n")
        f.write(f"ZONE T=\"Zone 1\", I={resolution}, J={resolution}, DATAPACKING=POINT\n")
        for j in range(resolution):
            for i in range(resolution):
                f.write(f"{i} {j} {u_mean[i, j]} {v_mean[i, j]} {u_std[i, j]} {v_std[i, j]}\n")


samples_number = 1000
total_time_cost = 137.54+135.969+143.776+142.977
#64 881.197+887.349
#128 8835.7
nx = ny =33
u0 = 0.01

u = np.zeros((samples_number, nx, ny))
v = np.zeros((samples_number, nx, ny))
tke = np.zeros(samples_number)

dir = "/home/zhongml95/github/sglbm/examples/tgv2d_mc/data/nx"+ str(nx)

tke_mean = np.zeros(samples_number)
tke_var = np.zeros(samples_number)

for n in range(0, samples_number):
    filename_u = dir+"/final/uAll_" + str(n) + ".dat"
    filename_v = dir+"/final/vAll_" + str(n) + ".dat"

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

print("tke_mean: ", tke_mean[-1])

u_mean = np.zeros((nx, ny))
v_mean = np.zeros((nx, ny))
u_std = np.zeros((nx, ny))
v_std = np.zeros((nx, ny))
tke_mean_all = 0

for i in range(0, nx):
    for j in range(0, ny):
        u_mean[i,j] = np.mean(u[:,i,j])
        v_mean[i,j] = np.mean(v[:,i,j])
        u_std[i,j] = np.std(u[:,i,j])
        v_std[i,j] = np.std(v[:,i,j])
        tke_mean_all += (u_mean[i,j]*u_mean[i,j] + v_mean[i,j]*v_mean[i,j]) * 2 / (nx*ny*u0*u0)


print("tke_mean_all: ", tke_mean_all)

save_to_tecplot(nx, u_mean, v_mean, u_std, v_std, dir + "/final/tgv2d.dat", title="Velocity Field")


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
plt.savefig("./data/computational_cost.png")



