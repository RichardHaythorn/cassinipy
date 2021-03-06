"""
File containing the plotting functions for CAPS data
"""
import heliopy.data.cassini as heliopydata
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


def base_els_spectrogram(starttime, endtime, anode, minmax, ax=None, cax=None, fig=None):
    """
    Plot spectrogram of ELS data

    Note: Not intended to be called by user

    :param anode:
    :param starttime:
    :param endtime:
    :param minmax:
    :param ax:
    :return:
    """
    capselsdata = heliopydata.caps_els(starttime, endtime, anode - 1, try_download=True)

    x = capselsdata.index
    y = range(capselsdata.shape[1])
    z = np.transpose(capselsdata.to_dataframe())

    maxvalue = capselsdata.to_dataframe().max().max()
    CS = ax.pcolormesh(x, y, z, norm=LogNorm(vmin=minmax[0], vmax=minmax[1]), cmap='jet')
    if fig is not None:
        cbar = fig.colorbar(cm.ScalarMappable(norm=LogNorm(vmin=minmax[0], vmax=minmax[1]), cmap='jet'), cax=cax)
        cbar.set_label("Counts/sec")

    ax.set_ylabel("CAPS ELS \nAnode {0} \nEnergy Bins ".format(anode))
    return CS, maxvalue


def els_spectrogram(starttime, endtime, anodes, minmax=(1, 2e6), subplot_kw={}):
    """
    Plot spectrogram of ELS data

    """
    if len(set(anodes).difference(set(range(1, 9, 1)))) is not 0:
        raise ValueError
    fig, axes = plt.subplots(nrows=len(anodes), subplot_kw=subplot_kw)
    specaxes = [axes] if not hasattr(axes, "__iter__") else axes
    for subplotcounter, anode in enumerate(anodes):
        base_els_spectrogram(starttime, endtime, anode, minmax, ax=specaxes[subplotcounter])

    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.9, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(cm.ScalarMappable(norm=LogNorm(vmin=minmax[0], vmax=minmax[1]), cmap='jet'), cax=cbar_ax)
    cbar.set_label("Counts/sec")
