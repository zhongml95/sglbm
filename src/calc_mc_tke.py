import numpy as np

resolution = 32
nMC = 200

tke = np.array(nMC)


for i in range(nMC):
    filenameU = "/home/zhongml95/sglbm/cluster/MC" + str(resolution) + "/u_" + str(i) + ".dat"
    u = np.loadtxt(filename)
    filenameU = "/home/zhongml95/sglbm/cluster/MC" + str(resolution) + "/u_" + str(i) + ".dat"
    v = np.loadtxt(filename)

    tke[i] = 0

    for nx in range(resolution)
        for ny in range(resolution)
            tke[i] = tke[i] + (u[nx][ny] * u[nx][ny] + v[nx][ny] * v[nx][ny])
    tke[i] = tke[i] * 2 / (nx * ny)