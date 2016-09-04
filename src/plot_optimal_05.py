import numpy as np
import matplotlib.pyplot as plt
import bicycleparameters as bp

# Plot Settings

figwidth_pt = 427.43153
figwidth = figwidth_pt / 72.27  # in inches
goldenMean = (np.sqrt(5.0) - 1.0) / 2.0
figsize = [figwidth, figwidth * goldenMean]
params = {'axes.labelsize': 8,
          'axes.titlesize': 8,
          'text.fontsize': 10,
          'legend.fontsize': 8,
          'xtick.labelsize': 6,
          'ytick.labelsize': 6,
          'figure.figsize': figsize,
          }
plt.rcParams.update(params)

speed = '05'
benchmark = bp.Bicycle('Benchmark',
                       pathToData='/home/moorepants/src/BicycleParameters/data')
optimal = bp.Bicycle('Optimal{}'.format(speed), pathToData='data')

fig = plt.figure()

axes = []
axes.append(plt.subplot2grid((2, 2), (0, 0), colspan=2))
axes.append(plt.subplot2grid((2, 2), (1, 0)))
axes.append(plt.subplot2grid((2, 2), (1, 1)))

# Geometry Plot

benchmark.plot_bicycle_geometry(show=False, inertiaEllipse=False, color='grey',
                                ax=axes[0])
optimal.plot_bicycle_geometry(show=False, inertiaEllipse=False, color='black',
                              ax=axes[0])
axes[0].set_title('Bicycle Geometry')
axes[0].set_xlabel('Distance [m]')
axes[0].set_ylabel('Distance [m]')

# Eigenvalues Versus Speed Plot

speeds = np.linspace(0.0, 10.0, num=100)

bp.plot_eigenvalues([benchmark, optimal], speeds, ax=axes[1],
                    colors=['grey', 'black'])

axes[1].set_ylim((-10, 15))
l = axes[1].legend(['Benchmark', 'Optimal @ 5 m/s'])
l.remove()
axes[1].set_ylabel('Real & Imaginary\nParts of the\nEigenvalues [1/s]')

# Handling Quality Metric Plot

benchmark_hqm_data = np.loadtxt('data/hqmbenchmark.csv'.format(speed),
                                delimiter=',')
optimal_hqm_data = np.loadtxt('data/hqm{}.csv'.format(speed), delimiter=',')

axes[2].plot(*benchmark_hqm_data.T, color='grey')
axes[2].plot(*optimal_hqm_data.T, color='black')
axes[2].set_title('Handling Quality Metric @ 5.0 m/s')
axes[2].set_xlabel('Frequency [rad/s]')
axes[2].set_ylabel('HQM')
axes[2].legend(['Benchmark', 'Optimal @ 5 m/s'])

plt.tight_layout()
plt.savefig('figures/example-optimal-bicycle.pdf')
