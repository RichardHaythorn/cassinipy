import heliopy.data.cassini as cassini
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


def base_ELS_spectrogram(anode, starttime, endtime, minmax, ax=None):
    """
    Plot spectrogram of ELS data
    :param anode:
    :param starttime:
    :param endtime:
    :param minmax:
    :param ax:
    :return:

    """
    capselsdata = cassini.caps_els(starttime, endtime, anode, try_download=True)

    x = capselsdata.index
    y = range(capselsdata.shape[1])
    z = np.transpose(capselsdata.to_dataframe())

    maxvalue = capselsdata.to_dataframe().max().max()
    CS = ax.pcolormesh(x, y, z, norm=LogNorm(vmin=minmax[0], vmax=minmax[1]), cmap='jet')
    return CS, maxvalue


def ELS_spectrogram(anodes, starttime, endtime, minmax=(1, 2e6), subplot_kw={}):
    """
    Plot spectrogram from SunPy Timeseries of ELS data

    """
    fig, specaxes = plt.subplots(nrows=len(anodes), subplot_kw=subplot_kw)
    for subplotcounter, anode in enumerate(anodes):
        base_ELS_spectrogram(anode, starttime, endtime, minmax, ax=specaxes[subplotcounter])
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(cm.ScalarMappable(norm=LogNorm(vmin=minmax[0], vmax=minmax[1]), cmap='jet'), cax=cbar_ax)
    cbar.set_label("Counts/sec")
