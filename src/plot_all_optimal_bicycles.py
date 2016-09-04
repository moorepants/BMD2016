import numpy as np
import matplotlib.pyplot as plt
import bicycleparameters as bp

# Plot Settings

figwidth_pt = 427.43153
figwidth = figwidth_pt / 72.27  # in inches
goldenMean = (np.sqrt(5.0) - 1.0) / 2.0
figsize = [figwidth, figwidth * 1.5]
params = {'axes.labelsize': 8,
          'axes.titlesize': 8,
          'text.fontsize': 10,
          'legend.fontsize': 8,
          'xtick.labelsize': 6,
          'ytick.labelsize': 6,
          'figure.figsize': figsize,
          }
plt.rcParams.update(params)

speed_strings = ['02', '04', '06', '08', '10']

fig = plt.figure()

hqm_ax = plt.subplot2grid((len(speed_strings), 2), (0, 1),
                          rowspan=len(speed_strings))
fig.subplots_adjust(hspace=0.0, wspace=0.0)

geom_axes = []

peak_hqm_values = []

for i, speed in enumerate(reversed(speed_strings)):

    bicycle = bp.Bicycle('Optimal{}'.format(speed), pathToData='data')

    if i > 0:
        geom_ax = plt.subplot2grid((len(speed_strings), 2), (i, 0),
                                sharex=geom_axes[-1])
        geom_ax.get_yaxis().set_visible(False)
    else:
        geom_ax = plt.subplot2grid((len(speed_strings), 2), (i, 0))

    bicycle.plot_bicycle_geometry(show=False, inertiaEllipse=False,
                                  ax=geom_ax)
                                  #centerOfMass=False, ax=geom_ax)

    geom_ax.set_ylim((0.0, 1.0))
    geom_ax.axis('off')
    geom_ax.set_title('')
    geom_axes.append(geom_ax)

    hqm_data = np.loadtxt('data/hqm{}.csv'.format(speed), delimiter=',')
    peak_hqm_values.append(hqm_data[:, 1].max())

peak_hqm_values.reverse()

hqm_ax.barh([float(s) - 1.0 for s in speed_strings],
            peak_hqm_values, height=2.0,
            log=True, color='grey')
hqm_ax.yaxis.tick_right()
hqm_ax.yaxis.set_label_position("right")
hqm_ax.set_ylim((1.0, 11.0))
hqm_ax.set_xlabel('max(HQM)')
hqm_ax.set_ylabel('Design Speed [m/s]')

plt.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.05)

fig.savefig('figures/all-optimal-bicycles.pdf')
