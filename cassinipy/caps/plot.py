import matplotlib.pyplot as plt
import numpy as np
from scipy.io import readsav
from sunpy.time import TimeRange

ibscalib = readsav('calib\\ibsdisplaycalib.dat')
elscalib = readsav('calib\\geometricfactor.dat')
sngcalib = readsav('calib\\sngdisplaycalib.dat')


def base_ELS_spectrogram(elstimeseries, starttime, endtime, ax=None, fig=None, subplot_kw={"yscale": "log"}):
    """
    Plot spectrogram from SunPy Timeseries of ELS data

    Returns:
    ax:
    """
    # TODO remove MSSL data dependency
    temptimeindex = elstimeseries.truncate(TimeRange(starttime, endtime)).index
    actualendtime = elstimeseries.data.index[elstimeseries.data.index.get_loc(temptimeindex[-1]) + 1]
    print(actualendtime)

    X = elstimeseries.truncate(TimeRange(starttime, actualendtime)).index
    Y = elscalib['polyearray']
    Z = np.transpose(elstimeseries.data.loc[starttime:endtime].iloc[:, 0:63].to_numpy())

    if ax is not None:
        specax = ax
        specfig = fig
    else:
        kwargs = subplot_kw
        specfig, specax = plt.subplots(subplot_kw=kwargs)

    # CS = specax.pcolormesh(X,Y,Z,norm=LogNorm(vmin=lower, vmax=upper),cmap='jet')
    CS = specax.pcolormesh(X, Y, Z, cmap='jet')
    cbar = specfig.colorbar(CS, ax=specax)
    # cbar.set_label("DEF [$m^{-2} s^{1} str^{-1} eV^{-1}$]")

    return specax


def ELS_spectrogram(elstimeseries, anodes, starttime, endtime, subplot_kw={"yscale": "log"}):
    """
    Plot spectrogram from SunPy Timeseries of ELS data

    Returns:
    ax:
    """
    specfig, specaxes = plt.subplots(nrows=len(anodes), subplot_kw=subplot_kw)
    for subplotcounter, anode in enumerate(anodes):
        base_ELS_spectrogram(elstimeseries[anode - 1], starttime, endtime, ax=specaxes[subplotcounter], fig=specfig)
