import numpy as np
import matplotlib.pyplot as plt




#diffusive scaling omega = 1.6
data8 = np.loadtxt("./t5/Nr4Nq9N8/final/tke.dat")
data16 = np.loadtxt("./t5/Nr4Nq9N16/final/tke.dat")
data32 = np.loadtxt("./t5/Nr4Nq9N32/final/tke.dat")
data64 = np.loadtxt("./t5/Nr4Nq9N64/final/tke.dat")
data128 = np.loadtxt("./t5/Nr5Nq11N128/final/tke.dat")

resolution = np.array([8,16,32,64])
err = np.zeros(resolution.size)
err[0] = abs(data8[0] - data128[0]) / data128[0]
err[1] = abs(data16[0] - data128[0]) / data128[0]
err[2] = abs(data32[0] - data128[0]) / data128[0]
err[3] = abs(data64[0] - data128[0]) / data128[0]

print(err)

# Fit a line to the data

coeffs = np.polyfit(np.log10(resolution), np.log10(err), deg=1)
fit_func = np.poly1d(coeffs)

fig, ax = plt.subplots()
#ax.plot(resolution, Cd, '-s', linewidth=2, markersize=8, markeredgewidth=1.5, markerfacecolor='none', label=r'$\bar{C_D}$')

#ax.loglog(resolution, err, linewidth=2, marker = "o", markersize=8, markeredgewidth=1.5, markerfacecolor='none', label='relative error, order={:.2f}'.format(coeffs[0]))

ax.scatter(resolution, err, marker = 'o', edgecolor = 'black',facecolor='none', label="error", s=30, zorder = 10, linewidths=2)

ax.loglog(resolution, 10**(fit_func(np.log10(resolution))), '-', 
           label='Fitted line, order={:.2f}'.format(coeffs[0]), linewidth = 2, color="#EA4335")

# Set axis labels and title
ax.set_xlabel('Resolution', fontsize=12)
ax.set_ylabel('relative error', fontsize=12)
#ax.set_title('relative error of drag coefficient at different resolutions')
ax.legend(loc='best', fontsize='large')

# Modify y-axis tick labels

# Show plot
ax.grid()
#plt.show()
fig.savefig('eoc.png')