import datetime

import matplotlib.pyplot as plt
import spiceypy as spice
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter
from numpy import array, max, arange, mean, transpose, ndarray, zeros, arccos
from pandas import DatetimeIndex, DataFrame
from scipy.io import readsav
from scipy.signal import find_peaks
from sunpy.time import TimeRange
from sunpy.timeseries import TimeSeriesMetaData, GenericTimeSeries

from cassinipy.caps.spice import cassini_ramdirection_SCframe, rotate_CAPS_SCframe

ibscalib = readsav('calib\\ibsdisplaycalib.dat')
elscalib = readsav('calib\\geometricfactor.dat')
sngcalib = readsav('calib\\sngdisplaycalib.dat')


def generate_timeseries_caps_mssl(data, anodefan, elsdatatype='data'):
    """
    Returns a SunPy Timeseries from processed MSSL data
    :param data:
    :param anodefan:
    :param elsdatatype:
    :return:
    """
    # TODO add multi fielded data?
    # TODO fix metadata

    if elsdatatype in data.keys():
        dataframe = DataFrame(transpose(data[elsdatatype][:, anodefan, :].byteswap().newbyteorder()))
        currentdate = datetime.datetime.strptime(data['sdate'].decode('ascii'), "%d-%b-%Y")
        for key, value in data.items():
            if isinstance(value, ndarray):
                if value.ndim == 1 and value.shape[0] == data[elsdatatype].shape[2]:
                    dataframe[key] = value.byteswap().newbyteorder()
        dataframe = dataframe.set_index(
            DatetimeIndex([currentdate + datetime.timedelta(seconds=x) for x in dataframe['secofday']]))
        tr = TimeRange(dataframe.index[0], dataframe.index[-1])
        dataframe.drop(columns=['time_ut', 'secofday', 'endsec', 'timehrs', 'hhmmss'], inplace=True)
        for i in range(63):
            dataframe.rename(columns={i: str(elscalib['earray'][i])}, inplace=True)

        caps_metadata = TimeSeriesMetaData(timerange=tr, colnames=['ELS'])
        # caps_metadata.append(tr,['ELS'],MetaDict([('formatid',dataframe['formatid'])]))

        capstimeseries = GenericTimeSeries(dataframe, meta=caps_metadata)

    # if 'sngdata' in data.keys():
    #     print("ims")
    #     instrument = "ims"
    # if 'ibsdata' in data.keys():
    #     print("ibs")
    #     instrument = "ibs"

    return capstimeseries


# data_metadata = TimeSeriesMetaData(meta=(instrument))

# return sunpy.timeseries.GenericTimeSeries(meta=data_metadata)


def CAPS_actuation(data, tempdatetime):
    """
    Returns the CAPS actuation give a datetime.datetime
    Works for elsres or sngres
    """

    for counter, i in enumerate(data['times_utc']):
        if i >= tempdatetime:
            slicenumber = counter
            break
    if 'sngact' in data.keys():
        actuation = data['sngact'][slicenumber]
    if 'actuator' in data.keys():
        dx = datetime.datetime.timestamp(data['times_utc'][slicenumber]) - datetime.datetime.timestamp(
            data['times_utc'][slicenumber - 1])
        dy = data['actuator'][slicenumber] - data['actuator'][slicenumber - 1]
        dy_dx = dy / dx
        newdx = datetime.datetime.timestamp(tempdatetime) - datetime.datetime.timestamp(
            data['times_utc'][slicenumber - 1])
        newdy = dy_dx * newdx
        actuation = data['actuator'][slicenumber - 1] + newdy

    return actuation


# TODO, make into function for ELS
# TODO, make into function for IBS

def CAPS_ELS_localramangle(tempdatetime, elsdata, anodes=False):
    act = CAPS_actuation(elsdata, tempdatetime)
    ramdir_SC = cassini_ramdirection_SCframe(tempdatetime, output=False)
    ELSvecs = rotate_CAPS_SCframe(act, 'els', anodes=anodes)

    if anodes:
        angle = zeros((8))
        for anodenumber, temp in enumerate(ELSvecs):
            angle[anodenumber] = arccos(spice.vdot(temp, ramdir_SC)) * spice.dpr()

    if not anodes:
        angle = arccos(spice.vdot(ELSvecs, ramdir_SC)) * spice.dpr()

    return angle


def CAPS_ELS_FOVcentre_azi_elv(tempdatetime, elsdata, anodes=False):
    """
    Either returns azi and elv for centre of anodes, or centre of full FOV
    """
    act = CAPS_actuation(elsdata, tempdatetime)
    ELSvecs = rotate_CAPS_SCframe(act, 'els', anodes=anodes)

    if anodes:
        ELS_azi = act
        ELS_elv = arange(-70, 90, 20)

    if not anodes:
        ELS_azi = act
        ELS_elv = 0

    return ELS_azi, ELS_elv


def CAPS_IBS_FOVcentre_azi_elv(tempdatetime, elsdata):
    """
    Returns azimuthal and elevation ram angles
    """
    # TODO : Remove dependence on ELS data
    # TODO: Add all fans

    act = CAPS_actuation(elsdata, tempdatetime)
    IBSvecs = rotate_CAPS_SCframe(act, 'ibs2', anodes=False)
    azimuthvec_norm = [0, -1, 0]
    elevationvec_norm = [0, 0, -1]

    # IBS_elv = np.arccos(spice.vdot(elevationvec_norm,IBSvecs))*spice.dpr()
    # IBS_azi = np.arcsin(spice.vdot(azimuthvec_norm,IBSvecs))*spice.dpr()
    ibs_elv = 0
    ibs_azi = act
    # print(act,IBSvecs,IBS_azi,IBS_elv)

    return ibs_azi, ibs_elv


def CAPS_actuationtimeslice(datetime, elsdata):
    """
    Returns the start and end slicenumbers while CAPS in actuation in one direction
    Accounts for extrema
    """
    # TODO remove ELSdata dependency
    slicevalue = CAPS_slicenumber(elsdata, datetime)

    if elsdata['actuator'][slicevalue + 2] > (elsdata['actuator'][slicevalue] + 0.5):
        direction = "positive"
    elif elsdata['actuator'][slicevalue + 2] < (elsdata['actuator'][slicevalue] - 0.5):
        direction = "negative"
    elif abs(elsdata['actuator'][slicevalue + 1] - elsdata['actuator'][slicevalue]) < 0.5:
        direction = "extrema"
    # print(direction)

    startactslice = slicevalue
    endactslice = slicevalue

    if direction == "positive":
        while elsdata['actuator'][startactslice - 2] < elsdata['actuator'][startactslice]:
            startactslice -= 1
        while elsdata['actuator'][endactslice + 2] > elsdata['actuator'][endactslice]:
            endactslice += 1

    if direction == "negative":
        while elsdata['actuator'][startactslice - 2] > elsdata['actuator'][startactslice]:
            startactslice -= 1
        while elsdata['actuator'][endactslice + 2] < elsdata['actuator'][endactslice]:
            endactslice += 1

    if direction == "extrema":
        while abs(elsdata['actuator'][startactslice - 1] - elsdata['actuator'][startactslice]) < 0.5:
            startactslice -= 1
        while abs(elsdata['actuator'][endactslice + 1] - elsdata['actuator'][endactslice]) < 0.5:
            endactslice += 1
        startactslice -= 10
        endactslice += 10

    return startactslice, endactslice, direction


def CAPS_slicenumber(data, tempdatetime):
    for counter, i in enumerate(data['times_utc']):
        if i >= tempdatetime:
            slicenumber = counter
            break

    return slicenumber


def CAPS_energyslice(sensor, startenergy, endenergy):
    """
    Returns the start and end slices of the earray required to cover energy range

    """

    if sensor == "els":
        polyearray = elscalib['polyearray']
    if sensor == "ims":
        polyearray = sngcalib['sngpolyearray']
    if sensor == "ibs":
        polyearray = ibscalib['ibspolyearray']

    if startenergy < polyearray[0]:
        startcounter = 0
    else:
        for counter, energy in enumerate(polyearray):
            if energy <= startenergy < polyearray[counter + 1]:
                startcounter = counter
                break

    if endenergy > polyearray[-1]:
        endcounter = len(polyearray) - 1
    else:
        for counter, energy in enumerate(polyearray):
            if energy <= endenergy < polyearray[counter + 1]:
                endcounter = counter
                break

    return startcounter, endcounter


def ELS_backgroundremoval(data, startslice, endslice):
    def_backgroundremoved = zeros((63, 8, endslice - startslice))
    # Background defined as average of 5 loweest count, negative ions unlikely to appear across 3 anodes
    for backgroundcounter, timecounter in enumerate(arange(startslice, endslice, 1)):

        for energycounter in range(63):
            backgroundremoved_temp = array(data['def'][energycounter, :8, timecounter]) - mean(
                sorted(data['def'][energycounter, :8, timecounter])[:5])
            backgroundremoved_anodes = [0 if i < 0 else i for i in backgroundremoved_temp]
            def_backgroundremoved[energycounter, :, backgroundcounter] = backgroundremoved_anodes

    return def_backgroundremoved


def ELS_intensitypeaks(elsdata, energy, anode, starttime, endtime, prominence=1e10):
    """
    Works on DEF
    :param elsdata:
    :param energy:
    :param anode:
    :param starttime:
    :param endtime:
    :param prominence:
    :return:
    """
    startslice, endslice = CAPS_slicenumber(elsdata, starttime), CAPS_slicenumber(elsdata, endtime)
    energybin = CAPS_energyslice("els", energy, energy)[0]

    dataslice = elsdata['def'][energybin, anode, startslice:endslice]
    peaks, properties = find_peaks(dataslice, height=1e3, prominence=prominence)
    times = []
    for x in peaks:
        times.append(elsdata['times_utc'][x + startslice])
    return times, properties


def IBS_intensitypeaks(ibsdata, energy, starttime, endtime, prominence=2e3):
    startslice, endslice = CAPS_slicenumber(ibsdata, starttime), CAPS_slicenumber(ibsdata, endtime)
    energybin = CAPS_energyslice("ibs", energy, energy)[0]

    dataslice = ibsdata['ibsdata'][energybin, 1, startslice:endslice]
    peaks, properties = find_peaks(dataslice, height=1e3, prominence=prominence)
    times = []
    for x in peaks:
        times.append(ibsdata['times_utc'][x + startslice])
    return times, properties


def ELS_intensityplot(elsdata, energy, anode, starttime, endtime, peaks=False, prominence=1e10,
                      subplot_kw={"yscale": "log"}):
    startslice, endslice = CAPS_slicenumber(elsdata, starttime), CAPS_slicenumber(elsdata, endtime)
    energybin = CAPS_energyslice("els", energy, energy)[0]

    fig, ax = plt.subplots(subplot_kw=subplot_kw)
    ax.plot(elsdata['times_utc'][startslice:endslice], elsdata['def'][energybin, anode - 1, startslice:endslice])
    ax.set_xlabel("Time")
    ax.set_ylabel("Counts")

    if peaks:
        times, properties = ELS_intensitypeaks(elsdata, energy, anode - 1, starttime, endtime, prominence=prominence)
        ax.scatter(times, properties['peak_heights'])


def IBS_intensityplot(ibsdata, energy, starttime, endtime, peaks=False, prominence=2e3, subplot_kw={"yscale": "log"}):
    startslice, endslice = CAPS_slicenumber(ibsdata, starttime), CAPS_slicenumber(ibsdata, endtime)
    energybin = CAPS_energyslice("ibs", energy, energy)[0]

    fig, ax = plt.subplots(subplot_kw=subplot_kw)
    ax.plot(ibsdata['times_utc'][startslice:endslice], ibsdata['ibsdata'][energybin, 1, startslice:endslice])
    ax.set_xlabel("Time")
    ax.set_ylabel("Counts")

    if peaks:
        times, properties = IBS_intensitypeaks(ibsdata, energy, starttime, endtime, prominence=prominence)
        ax.scatter(times, properties['peak_heights'])


def ELS_spectrogram(elsdata, anodes, starttime, seconds, maxbackgroundremoved=False, meanbackgroundremoved=False):
    """
    Plot multiple els anodes from single bin across time
    """

    deflimits = {'t27': [1e9, 1e12], 't55': [1e9, 1e12], 't56': [1e9, 1e12], 't57': [1e9, 1e12], 't58': [1e9, 1e12],
                 't59': [1e9, 1e12]}

    energylimits = {'t27': [2, 500], 't40': [1, 10000], 't55': [2, 500], 't56': [2, 500], 't57': [2, 2e3],
                    't58': [2, 500], 't59': [2, 500], 't83': [2, 300]}

    if len(anodes) == 1:
        fig = plt.figure(figsize=(18, 6))
        axis1 = plt.axes()

    if len(anodes) > 1:
        fig, axes = plt.subplots(len(anodes), sharex=True, sharey=True, figsize=(18, 10))

    for counter, i in enumerate(elsdata['times_utc_strings']):
        if i >= starttime:
            slicenumber = counter
            actualtime = i
            break

    print(slicenumber)

    slicenumbers = arange(slicenumber, slicenumber + (seconds / 2), 1, dtype=int)

    X = elsdata['times_utc'][slicenumbers[0]:slicenumbers[-1] + 1]
    # print(X)
    Y = elscalib['earray']

    for figcounter, anode in enumerate(anodes):

        if maxbackgroundremoved and meanbackgroundremoved:
            Z = elsdata['def'][:, anode, slicenumbers[0]:slicenumbers[-1]]
        if meanbackgroundremoved:
            backgroundremoveddata = zeros(shape=(63, len(X)))
            for slicecounter, slicenumber in enumerate(slicenumbers):
                backgroundcount = mean([mean(elsdata['def'][:, :3, slicenumber], axis=1),
                                        mean(elsdata['def'][:, 5:8, slicenumber], axis=1)], axis=0)
                tempdata = [a - b for a, b in zip(elsdata['def'][:, anode, slicenumber], backgroundcount)]
                backgroundremoveddata[:, slicecounter] = tempdata
            Z = backgroundremoveddata
        if maxbackgroundremoved:
            backgroundremoveddata = zeros(shape=(63, len(X)))
            for slicecounter, slicenumber in enumerate(slicenumbers):
                backgroundcount = max([max(elsdata['def'][:, :3, slicenumber], axis=1),
                                       max(elsdata['def'][:, 5:8, slicenumber], axis=1)], axis=0)
                tempdata = [a - b for a, b in zip(elsdata['def'][:, anode, slicenumber], backgroundcount)]
                backgroundremoveddata[:, slicecounter] = tempdata
            Z = backgroundremoveddata

        print(Z.shape)

        if elsdata['flyby'] in deflimits.keys():
            lower = deflimits[elsdata['flyby']][0]
            upper = deflimits[elsdata['flyby']][1]
        else:
            lower = 1e8
            upper = 1e14
        if len(anodes) > 1:
            axis1 = axes[figcounter]
        cs = axis1.pcolormesh(X, Y, Z, norm=LogNorm(vmin=lower, vmax=upper), cmap='jet')

        axis1.set_yscale("log")
        axis1.set_xlim(X[0], X[-1])
        if elsdata['flyby'] in energylimits.keys():
            axis1.set_ylim(energylimits[elsdata['flyby']][0], energylimits[elsdata['flyby']][1])
        axis1.set_ylabel("Anode " + str(anode + 1) + "\n [eV/q]", fontsize=20)
        axis1.minorticks_on()
        axis1.tick_params(labelbottom=True, labeltop=False, bottom=True, top=True, left=True, right=True, which='both')
        axis1.yaxis.set_major_formatter(ScalarFormatter())

        if len(anodes) < 4:
            axis1.tick_params(axis='y', labelleft=True, labelright=False, left=True, right=True, which='both')

        if len(anodes) >= 4:
            axis1.tick_params(axis='y', labelleft=True, labelright=False, left=True, right=True, which='both',
                              labelsize=10)

            if anode != anodes[-1]:
                axis1.set_xticklabels([''])

    if len(anodes) == 1:
        axis1.set_title("Spectrogram of " + elsdata['instrument'].upper() + " data from the " + str(
            elsdata['flyby']).upper() + " flyby", y=1.01, fontsize=32)
        axis1.set_xlabel("Time", fontsize=20)
    if len(anodes) > 1:
        axes[0].set_title("Spectrogram of " + elsdata['instrument'].upper() + " data from the " + str(
            elsdata['flyby']).upper() + " flyby", y=1.01, fontsize=32)
        axes[-1].set_xlabel("Time", fontsize=20)

    fig.autofmt_xdate()

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.82, 0.2, 0.03, 0.68])
    cbar = fig.colorbar(cs, cax=cbar_ax)
    cbar.ax.set_ylabel("DEF [$m^{-2} s^{1} str^{-1} eV^{-1}$]")
    # fig.savefig(elsdata['instrument'] + "-" + elsdata['flyby'] + "-Spectrogram.pdf",format='pdf',bbox_inches='tight')

    plt.show()


def SNG_spectrogram(sngdata, anodes, starttime, seconds, mass=False, backgroundremoved=False, hlines=False,
                    hlines_vel=False):
    """
    Plot multiple sng anodes from single bin across time
    """
    deflimits = {'e3': [1e10, 1e12], 'e5': [1e11, 1e14], 'e7': [1e9, 1e13], 'e17': [1e9, 5e12], 'e18': [1e9, 5e12]}

    energylimits = {'e5': [10, 70], 'e7': [1, 70], 'e17': [1, 70], 'e18': [1, 70]}
    masslimits = {'e5': [0.5, 0.45e2], 'e7': [10, 0.5e2], 'e17': [10, 200]}

    if len(anodes) == 1:
        fig = plt.figure(figsize=(18, 6))
        axis1 = plt.axes()

    if len(anodes) > 1:
        fig, axes = plt.subplots(len(anodes), sharex='True', sharey='True', figsize=(18, 10))

    for counter, i in enumerate(sngdata['times_utc_strings']):
        if i >= starttime:
            slicenumber = counter
            actualtime = i
            break

    slicenumbers = arange(slicenumber, slicenumber + (seconds / 4), 1, dtype=int)
    print(slicenumbers, len(slicenumbers))
    print(slicenumbers[0], slicenumbers[-1])

    X = sngdata['times_utc'][slicenumbers[0]:slicenumbers[-1] + 1]

    if not mass:
        Y = sngcalib['sngpolyearray']
    if mass:
        averageSCP = -6.5
        print(averageSCP)
        Y = [(i - averageSCP) * sngdata['conversionfactor'] for i in sngcalib['sngpolyearray']]

    for figcounter, anode in enumerate(anodes):

        if not backgroundremoved:
            Z = sngdata['sngdef'][:, anode, slicenumbers[0]:slicenumbers[-1]]
        if backgroundremoved:
            backgroundremoveddata = zeros(shape=(62, len(slicenumbers)))
            for slicecounter, slicenumber in enumerate(slicenumbers):
                backgroundcount = mean([mean(sngdata['sngdef'][:, :3, slicenumber], axis=1),
                                        mean(sngdata['sngdef'][:, 5:8, slicenumber], axis=1)], axis=0)
                tempdata = [a - b for a, b in zip(sngdata['sngdef'][:, anode, slicenumber], backgroundcount)]
                backgroundremoveddata[:, slicecounter] = tempdata
            Z = backgroundremoveddata

        print(len(X), len(Y), Z.shape)

        if len(anodes) > 1:
            axis1 = axes[figcounter]
        if 'flyby' in list(deflimits.keys()):
            cs = axis1.pcolormesh(X, Y, Z, norm=LogNorm(vmin=deflimits[sngdata['flyby']][0],
                                                        vmax=deflimits[sngdata['flyby']][1]), cmap='jet')
            axis1.set_ylim(energylimits[sngdata['flyby']][0], energylimits[sngdata['flyby']][1])
        else:
            cs = axis1.pcolormesh(X, Y, Z, norm=LogNorm(vmin=5e10, vmax=5e14), cmap='jet')

        axis1.set_yscale("log")
        axis1.set_xlim(X[0], X[-1])

        axis1.set_ylabel("A" + str(anode + 1) + "\n [/eV]", fontsize=20)
        axis1.minorticks_on()

    if len(anodes) == 1:
        axis1.set_title("Spectrogram of IMS data from the " + str(sngdata['flyby']).upper() + " flyby", y=1.01,
                        fontsize=32)
        axis1.set_xlabel("Time", fontsize=20)
    if len(anodes) > 1:
        axes[0].set_title("Spectrogram of IMS data from the " + str(sngdata['flyby']).upper() + " flyby", y=1.01,
                          fontsize=32)
        axes[-1].set_xlabel("Time", fontsize=20)

    fig.autofmt_xdate()

    if mass == True or hlines == True or hlines_vel == True:
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.89, 0.2, 0.02, 0.68])
        cbar = fig.colorbar(cs, cax=cbar_ax)
        cbar.ax.set_ylabel("DEF [$m^{-2} s^{1} str^{-1} eV^{-1}$]")
    else:
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.82, 0.2, 0.03, 0.68])
        cbar = fig.colorbar(cs, cax=cbar_ax)
        cbar.ax.set_ylabel("DEF [$m^{-2} s^{1} str^{-1} eV^{-1}$]")
    # fig.savefig(sngdata['instrument'] + "-" + sngdata['flyby'] + "-IBSSpectrogram.pdf",format='pdf',bbox_inches='tight')


def IBS_spectrogram(ibsdata, fans, starttime, seconds, backgroundremoved=False):
    """
    Plot multiple ibs fan spectrograms
    """

    countlimits = {'t27': [1e1, 5e5], 't46': [1e1, 5e5], 't55': [1e1, 5e5], 't56': [1e1, 5e5], 't57': [1e2, 5e5],
                   't58': [1e1, 5e5], 't59': [6e2, 3e5]}
    deflimits = {'t27': [1e10, 1e14], 't46': [1e10, 1e14], 't55': [1e10, 1e14], 't56': [1e10, 1e14],
                 't57': [1e10, 1e14], 't58': [1e10, 1e14], 't59': [1e10, 1e14]}

    energylimits = {'t27': [3, 100], 't40': [1, 10000], 't46': [3, 100], 't55': [35, 60], 't56': [35, 70],
                    't57': [1, 80], 't58': [35, 70], 't59': [3, 60], 't83': [1, 250]}

    if len(fans) == 1:
        fig = plt.figure(figsize=(18, 6))
        axis1 = plt.axes()

    if len(fans) > 1:
        fig, axes = plt.subplots(len(fans), sharex="True", sharey="True", figsize=(18, 10))

    for counter, i in enumerate(ibsdata['times_utc_strings']):
        if i >= starttime:
            slicenumber = counter
            break

    slicenumbers = arange(slicenumber, slicenumber + (seconds / 2), 1, dtype=int)

    X = ibsdata['times_utc'][slicenumbers[0]:slicenumbers[-1] + 1]
    Y = ibscalib['ibspolyearray']

    for figcounter, fan in enumerate(fans):

        if not backgroundremoved:
            Z = zeros((653, len(slicenumbers)))
            for k in range(653):
                for slicecounter, slicenumber in enumerate(slicenumbers):
                    Z[k, slicecounter] = ((ibsdata['ibsdata'][k, fan, slicenumber]) / (ibscalib['ibsgeom'] * 1e-4))
                    # Z[k,slicecounter] = ibsdata['ibsdata'][k,fan,slicenumber]
        if backgroundremoved:
            backgroundremoveddata = zeros(shape=(653, len(slicenumbers)))
            for slicecounter, slicenumber in enumerate(slicenumbers):
                backgroundval = max([max(ibsdata['ibsdata'][:, :3, slicenumber], axis=1),
                                     max(ibsdata['ibsdata'][:, 5:8, slicenumber], axis=1)], axis=0)
                tempdata = [a - b for a, b in zip(ibsdata['ibsdata'][:, fan, slicenumber], backgroundcount)]
                backgroundremoveddata[:, slicecounter] = tempdata
            Z = backgroundremoveddata

        print(len(X), len(Y), Z.shape)
        if len(fans) > 1:
            axis1 = axes[figcounter]

        if 'flyby' in list(deflimits.keys()):
            CS = axis1.pcolormesh(X, Y, Z, norm=LogNorm(vmin=deflimits[ibsdata['flyby']][0],
                                                        vmax=deflimits[ibsdata['flyby']][1]), cmap='jet')
            axis1.set_ylim(energylimits[ibsdata['flyby']][0], energylimits[ibsdata['flyby']][1])
        else:
            CS = axis1.pcolormesh(X, Y, Z, norm=LogNorm(vmin=5e10, vmax=5e14), cmap='jet')

        axis1.set_yscale("log")
        axis1.set_ylabel("IBS Fan " + str(fan + 1) + "\n [eV/q]", fontsize=16)
        axis1.minorticks_on()
        axis1.tick_params(labelbottom=True, labeltop=False, bottom=True, top=True, left=True, right=True, which='both')
        axis1.tick_params(axis='y', labelleft=True, labelright=False, left=True, right=True, which='major')
        axis1.yaxis.set_major_formatter(ScalarFormatter())
        # axis1.yaxis.set_minor_formatter(ScalarFormatter())

    if 'flyby' in list(ibsdata.keys()):
        if len(fans) == 1:
            axis1.set_title("Spectrogram of " + ibsdata['instrument'].upper() + " data from the " + str(
                ibsdata['flyby']).upper() + " flyby", y=1.01, fontsize=30)
            axis1.set_xlabel("Time", fontsize=20)
        if len(fans) > 1:
            axes[0].set_title("Spectrogram of " + ibsdata['instrument'].upper() + " data from the " + str(
                ibsdata['flyby']).upper() + " flyby", y=1.08, fontsize=30)
            axes[-1].set_xlabel("Time", fontsize=20)

    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.91, 0.11, 0.03, 0.77])
    cbar = fig.colorbar(CS, cax=cbar_ax)
    cbar.ax.set_ylabel("DEF [$m^{-2} s^{1} str^{-1} eV^{-1}$]")
    # cbar.ax.set_ylabel("Counts [/s]")
    # fig.savefig(ibsdata['instrument'] + "-" + ibsdata['flyby'] + "-IBSSpectrogram.pdf",format='pdf',bbox_inches='tight')
