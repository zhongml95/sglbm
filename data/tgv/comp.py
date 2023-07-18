import numpy as np
import matplotlib.pyplot as plt



SGdataU = np.loadtxt("./t5/Nr4Nq9N32/U1646.dat")

mcDataU1E3 = np.loadtxt("/home/zhongml95/olb/examples/laminar/tgv2d/results/NX32/uniform/N1/Nq100/tgv2dU.dat")
mcDataU1E4 = np.loadtxt("/home/zhongml95/olb/examples/laminar/tgv2d/results/NX32/uniform/N1/Nq1000/tgv2dU.dat")
mcDataU1E5 = np.loadtxt("/home/zhongml95/olb/examples/laminar/tgv2d/results/NX32/uniform/N1/Nq10000/tgv2dU.dat")
#mcDataU1E6 = np.loadtxt("/home/zhongml95/olb/examples/laminar/tgv2d/results/NX32/uniform/N1/Nq100000/tgv2dU.dat")

colors = ['#EA4335', '#4285F4', '#FBBC05', '#34A853', '#7B0099', '#00C853']

# umean
plt.figure(1,dpi=400)
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
plt.plot(mcDataU1E3[:,0], mcDataU1E3[:,1], label="MCS(OpenLB, 1E3)", linewidth = 3, color=colors[1], linestyle='-')
plt.plot(SGdataU[:,0], SGdataU[:,1], label="SG", linewidth = 3, color=colors[0],linestyle='--')
#plt.plot(SGdataU[:,0], SGdataU[:,3], label="Exact", linewidth = 3, color=colors, linestyle=':')
plt.scatter(SGdataU[:,0], SGdataU[:,3], marker = 'o', edgecolor = 'black',facecolor='none', label="Exact", s=30, zorder = 10, linewidths=2)
plt.rcParams.update({'font.size': 12})
# Setting x-axis as log scale
plt.legend(fontsize = 'large')
plt.grid(True)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# Setting labels and title
plt.xlabel('Y', fontsize=14)
plt.ylabel(r'$\bar{u}$', fontsize=14, labelpad=0)

# Display the plot
plt.savefig("./tgvMean.png")

# uStd
plt.figure(2,dpi=400)
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
#plt.plot(SGdataU[:,0], SGdataU[:,2], label="SG", linewidth = 3, color=colors[0],linestyle='-')
plt.plot(mcDataU1E3[:,0], mcDataU1E3[:,2], label="MCS(OpenLB, 1E3)", linewidth = 2, color=colors[1], linestyle='-')
plt.plot(mcDataU1E4[:,0], mcDataU1E4[:,2], label="MCS(OpenLB, 1E4)", linewidth = 2, color=colors[2], linestyle='-')
plt.plot(mcDataU1E5[:,0], mcDataU1E5[:,2], label="MCS(OpenLB, 1E5)", linewidth = 2, color=colors[3], linestyle='-')
#plt.plot(mcDataU1E6[:,0], mcDataU1E6[:,2], label="MC(OpenLB, 1E6)", linewidth = 2, color=colors[4], linestyle='-')
plt.plot(SGdataU[:,0], SGdataU[:,2], label="SG", linewidth = 2, color=colors[0],linestyle='--')
plt.rcParams.update({'font.size': 12})
# Setting x-axis as log scale
plt.legend(fontsize = 'large')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

plt.grid(True)

# Setting labels and title
plt.xlabel('Y', fontsize=14)
plt.ylabel(r'${\sigma}_{u}$', fontsize=14, labelpad=0)

# Display the plot
plt.savefig("./tgvStd.png")
