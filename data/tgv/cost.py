import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


def log_tick_formatter(val, pos=None):
    return f"$10^{{{int(val)}}}$"

Nq = np.array([10, 20, 40, 80, 160, 320, 640])

timeNr1 = np.array([7.03431, 8.15733, 10.7976, 15.5964, 22.8505, 41.7675, 72.5611])
timeNr2 = np.array([8.00459, 8.67968, 12.0963, 20.0002, 34.7041, 57.4981, 95.907])
timeNr3 = np.array([9.94482, 11.4869, 15.3234, 22.6927, 38.7666, 66.4299, 116.59])
timeNr4 = np.array([12.0751, 14.9576, 17.4408, 25.474, 47.6208, 76.3714, 141.677])
timeNr5 = np.array([14.3943, 17.0025, 23.19  , 30.6632, 51.5601, 89.7594, 167.539])


fig = plt.figure(dpi=400)
ax = plt.axes()

#ax.zaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
#ax.zaxis.set_major_locator(mticker.MaxNLocator(integer=True))

#ax.set_xscale('log', base=2)

ax.plot(Nq, timeNr1, marker='o', label=r'$N=1$')
ax.plot(Nq, timeNr2, marker='o', label=r'$N=2$')
ax.plot(Nq, timeNr3, marker='o', label=r'$N=3$')
ax.plot(Nq, timeNr4, marker='o', label=r'$N=4$')
ax.plot(Nq, timeNr5, marker='o', label=r'$N=5$')

#ax.legend(ncols=4, loc=8, bbox_to_anchor=(0.5, -0.15), fontsize='small')

ax.set_xlabel(r'Quadrature points number $M$', fontsize=12)
plt.legend(fontsize = 'large')

#ax.set_ylabel('Polynomial Order')
ax.set_ylabel('time cost', fontsize=12)

fig.savefig("/home/zhongml95/github/sglbm/data/tgv/sgCostNq.png")

fig = plt.figure(dpi=400)
ax = plt.axes()

#ax.zaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
#ax.zaxis.set_major_locator(mticker.MaxNLocator(integer=True))

#ax.set_xscale('log', base=2)

ax.plot(np.array([1,2,3,4,5]), np.array([72.5611, 95.907, 116.59, 141.677, 167.539]), marker='o', label=r'$M=640$')

#ax.legend(ncols=4, loc=8, bbox_to_anchor=(0.5, -0.15), fontsize='small')

plt.legend(fontsize = 'large')

ax.set_xlabel('Polynomial Order n', fontsize=12)
#ax.set_ylabel('Polynomial Order')
ax.set_ylabel('time cost', fontsize=12)

fig.savefig("/home/zhongml95/github/sglbm/data/tgv/sgCostNr.png")