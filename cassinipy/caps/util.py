from datetime import datetime

import h5py as h5py
import numpy as np
import pandas as pd
from scipy import constants
from sunpy.timeseries import GenericTimeSeries

AMU = constants.physical_constants['atomic mass constant'][0]
AMU_eV = constants.physical_constants['atomic mass unit-electron volt relationship'][0]
e = constants.physical_constants['atomic unit of charge'][0]
e_mass = constants.physical_constants['electron mass'][0]
e_mass_eV = constants.physical_constants['electron mass energy equivalent in MeV'][0] * 1e6
c = constants.physical_constants['speed of light in vacuum'][0]
k = constants.physical_constants['Boltzmann constant'][0]


def caps_dateparser(stringlist):
    return [datetime.strptime(item, "%Y-%jT%H:%M:%S.%f") for item in stringlist]


def caps_hdf5totimeseries(hdfpath, instrumentanode=0):
    """
    Creates a timeseries for a single anode/fan for a CAPS sensor
    """
    # TODO add the metadata
    # TODO add the units
    hdf5file = h5py.File(hdfpath, 'r')
    data_df = pd.DataFrame(hdf5file['DATA'][:, :, instrumentanode])
    times = caps_dateparser([x.decode("ASCII") for x in hdf5file['UTC']])
    data_df = data_df.set_index(pd.DatetimeIndex(times))
    date_timeseries = GenericTimeSeries(data_df)
    hdf5file.close()
    return date_timeseries


def peakflux(mass, spacecraftvelocity, flowspeed, spacecraftpotential, temperature, deflectionvelocity=0,
             electricfield=0, charge=1):
    """
    (mass,spacecraftvelocity,flowspeed,spacecraftpotential,temperature,charge)
    """

    mass_kg = mass * AMU

    spacecraftsurface = 2

    if electricfield != 0:
        deflectionvelocity = (e * electricfield * spacecraftsurface) / (mass_kg * spacecraftvelocity)
        print(deflectionvelocity)

    peakenergy = 0.5 * mass_kg * ((
                                          spacecraftvelocity + flowspeed + charge * deflectionvelocity) ** 2) - charge * spacecraftpotential * e + 8 * k * temperature
    peakenergy_ev = peakenergy / e

    return peakenergy_ev


def energy2mass(energy, spacecraftvelocity, flowspeed, spacecraftpotential, temperature, charge=1):
    """
    Energy in eV
    spacecraftpotential in volts
    temperature in K
    spacecraftvelocity in m/s
    flowspeed in m/s
    """
    # print(energy,spacecraftpotential,temperature,spacecraftvelocity,flowspeed)

    mass = (((energy * e + charge * spacecraftpotential * e - 8 * k * temperature) * 2) / (
            (spacecraftvelocity + flowspeed) ** 2)) / AMU

    return mass


def thermalenergy(temperature):
    thermalenergy_ev = (8 * k * temperature) / e

    return thermalenergy_ev


def thermalvelocity(mass, temperature):
    velocity = np.sqrt((k * temperature) / (mass * AMU))

    return velocity
