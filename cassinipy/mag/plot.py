"""
File for generating MAG data plots
"""
import heliopy.data.cassini as heliopydata


def base_mag_1min_plot(starttime, endtime, coords, ax=None):
    """
    Plot spectrogram of MAG data
    """
    magdata = heliopydata.mag_1min(starttime, endtime, coords)
    for Baxis in ['Bx', 'By', 'Bz']:
        ax.plot(magdata.index, magdata.to_dataframe()[Baxis], label=Baxis)
    ax.set_ylabel("MAG \n{0} coords  \n[\\nT] ".format(coords))
