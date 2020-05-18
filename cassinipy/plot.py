"""
Contains the main plotting functions for cassinipy
"""
from cassinipy.caps.plot import *
from cassinipy.mag.plot import *


def cassini_plot(starttime, endtime, **kwargs):
    """
    Main plotting function for cassinipy

    Keyword Arguments:
    caps_els = list of wanted anodes
    mag_1min = coords system

    :param starttime:
    :param endtime:
    :param subplot_kw:
    :param kwargs:
    """
    numberofsubplots = 0
    for value in kwargs.values():
        if isinstance(value, list):
            numberofsubplots += len(value)
        if isinstance(value, str):
            numberofsubplots += 1

    fig = plt.figure(constrained_layout=True)
    spec = fig.add_gridspec(ncols=2, nrows=numberofsubplots, width_ratios=[20, 1])

    subplotcounter = 0
    for key, value in kwargs.items():
        if key == "caps_els":
            if len(set(value).difference(set(range(1, 9, 1)))) is not 0:
                raise ValueError("Invalid anode number given, must be in {}".format(set(range(1, 9, 1))))
            for anode in value:
                ax = fig.add_subplot(spec[subplotcounter, 0])
                cax = fig.add_subplot(spec[subplotcounter, 1])
                base_els_spectrogram(starttime, endtime, anode, minmax=(1, 2e6), ax=ax, cax=cax,
                                     fig=fig)
                subplotcounter += 1
        if key == "mag_1min":
            ax = fig.add_subplot(spec[subplotcounter, 0])
            base_mag_1min_plot(starttime, endtime, value, ax=ax)
            subplotcounter += 1
