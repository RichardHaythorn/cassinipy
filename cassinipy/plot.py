"""
Contains the main plotting functions for cassinipy
"""
from cassinipy.caps.plot import *


def cassini_plot(starttime, endtime, subplot_kw={}, **kwargs):
    """

    :param starttime:
    :param endtime:
    :param subplot_kw:
    :param kwargs:
    """
    numberofsubplots = 0
    for value in kwargs.values():
        numberofsubplots += len(value)
    fig, tempaxes = plt.subplots(nrows=numberofsubplots, subplot_kw=subplot_kw)
    axes = [tempaxes] if not hasattr(tempaxes, "__iter__") else tempaxes
    subplotcounter = 0
    for key, value in kwargs.items():
        if key == "caps_els":
            if len(set(value).difference(set(range(1, 9, 1)))) is not 0:
                raise ValueError
            for anode in value:
                base_ELS_spectrogram(starttime, endtime, anode, minmax=(1, 2e6), ax=axes[subplotcounter],
                                     fig=fig)
                subplotcounter += 1
