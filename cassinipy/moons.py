from datetime import datetime
import spiceypy as spice
from collections import OrderedDict

flybydict = {'T0': ['titan', '2004-07-03'], 'TA': ['titan', '2004-10-26'], 'TB': ['titan', '2004-12-13'],
             'TC': ['titan', '2005-01-13'], 'T3': ['titan', '2005-02-15'], 'T4': ['titan', '2005-03-31'],
             'T5': ['titan', '2005-04-16'], 'T6': ['titan', '2005-08-22'], 'T7': ['titan', '2005-09-07'],
             'T8': ['titan', '2005-10-28'], 'T9': ['titan', '2005-12-26'], 'T10': ['titan', '2006-01-15'],
             'T11': ['titan', '2006-02-27'], 'T12': ['titan', '2006-03-19'], 'T13': ['titan', '2006-04-30'],
             'T14': ['titan', '2006-05-20'], 'T15': ['titan', '2006-07-02'], 'T16': ['titan', '2006-07-22'],
             'T17': ['titan', '2006-09-07'], 'T18': ['titan', '2006-09-23'], 'T19': ['titan', '2006-10-09'],
             'T20': ['titan', '2006-10-25'], 'T21': ['titan', '2006-12-12'], 'T22': ['titan', '2006-12-28'],
             'T23': ['titan', '2007-01-13'], 'T24': ['titan', '2007-01-29'], 'T25': ['titan', '2007-02-22'],
             'T26': ['titan', '2007-03-10'], 'T27': ['titan', '2007-03-26'], 'T28': ['titan', '2007-04-10'],
             'T29': ['titan', '2007-04-26'], 'T30': ['titan', '2007-05-12'], 'T31': ['titan', '2007-05-28'],
             'T32': ['titan', '2007-06-13'], 'T33': ['titan', '2007-06-29'], 'T34': ['titan', '2007-07-19'],
             'T35': ['titan', '2007-08-31'], 'T36': ['titan', '2007-10-02'], 'T37': ['titan', '2007-11-19'],
             'T38': ['titan', '2007-12-05'], 'T39': ['titan', '2007-12-20'], 'T40': ['titan', '2008-01-05'],
             'T41': ['titan', '2008-02-22'], 'E3': ['enceladus', '2008-03-12'], 'T42': ['titan', '2008-03-25'],
             'T43': ['titan', '2008-05-12'], 'T44': ['titan', '2008-05-28'], 'T45': ['titan', '2008-07-31'],
             'T46': ['titan', '2008-11-03'], 'T47': ['titan', '2008-11-19'], 'T48': ['titan', '2008-12-05'],
             'T49': ['titan', '2008-12-21'], 'T50': ['titan', '2009-02-07'], 'T51': ['titan', '2009-03-27'],
             'T52': ['titan', '2009-04-04'], 'T53': ['titan', '2009-04-20'], 'T54': ['titan', '2009-05-05'],
             'T55': ['titan', '2009-05-21'], 'T56': ['titan', '2009-06-06'], 'T57': ['titan', '2009-06-22'],
             'T58': ['titan', '2009-07-08'], 'T59': ['titan', '2009-07-24'], 'T60': ['titan', '2009-08-09'],
             'T61': ['titan', '2009-08-25'], 'T62': ['titan', '2009-10-12'], 'T63': ['titan', '2009-12-12'],
             'T64': ['titan', '2009-12-28'], 'T65': ['titan', '2010-01-12'], 'T66': ['titan', '2010-01-28'],
             'T67': ['titan', '2010-04-05'], 'T68': ['titan', '2010-05-20'], 'T69': ['titan', '2010-06-05'],
             'T70': ['titan', '2010-06-21'], 'T71': ['titan', '2010-07-07'], 'T72': ['titan', '2010-09-24'],
             'T73': ['titan', '2010-11-11'], 'T74': ['titan', '2011-02-18'], 'T75': ['titan', '2011-04-19'],
             'T76': ['titan', '2011-05-08'], 'T77': ['titan', '2011-06-20'], 'T78': ['titan', '2011-09-12'],
             'T79': ['titan', '2011-12-13'], 'T80': ['titan', '2012-01-02'], 'T81': ['titan', '2012-01-30'],
             'T82': ['titan', '2012-02-19'], 'T83': ['titan', '2012-05-22'], 'T84': ['titan', '2012-06-06'],
             'T85': ['titan', '2012-07-24'], 'T86': ['titan', '2012-09-26'], 'T87': ['titan', '2012-11-13'],
             'T88': ['titan', '2012-11-29'], 'T89': ['titan', '2013-02-17'], 'T90': ['titan', '2013-04-05'],
             'T91': ['titan', '2013-05-23'], 'T92': ['titan', '2013-07-10'], 'T93': ['titan', '2013-07-26'],
             'T94': ['titan', '2013-09-12'], 'T95': ['titan', '2013-10-14'], 'T96': ['titan', '2013-12-01'],
             'T97': ['titan', '2014-01-01'], 'T98': ['titan', '2014-02-02'], 'T99': ['titan', '2014-03-06'],
             'T100': ['titan', '2014-04-07'], 'T101': ['titan', '2014-05-17'], 'T102': ['titan', '2014-06-18'],
             'T103': ['titan', '2014-07-20'], 'T104': ['titan', '2014-08-21'], 'T105': ['titan', '2014-09-22'],
             'T106': ['titan', '2014-10-24'], 'T107': ['titan', '2014-12-10'], 'T108': ['titan', '2015-01-11'],
             'T109': ['titan', '2015-02-12'], 'T110': ['titan', '2015-03-16'], 'T111': ['titan', '2015-05-07'],
             'T112': ['titan', '2015-07-07'], 'T113': ['titan', '2015-09-28'], 'T114': ['titan', '2015-11-13'],
             'T115': ['titan', '2016-01-16'], 'T116': ['titan', '2016-02-01'], 'T117': ['titan', '2016-02-16'],
             'T118': ['titan', '2016-04-04'], 'T119': ['titan', '2016-05-06'], 'T120': ['titan', '2016-06-07'],
             'T121': ['titan', '2016-07-25'], 'T122': ['titan', '2016-08-10'], 'T123': ['titan', '2016-09-27'],
             'T124': ['titan', '2016-11-14'], 'T125': ['titan', '2016-11-29'], 'T126': ['titan', '2017-04-22']}


class MoonFlyby:
    def __init__(self, flyby, instruments):
        self.flyby = str(flyby)
        self.moon = flybydict[str(flyby)][0]
        self.date = datetime.strptime(flybydict[str(flyby)][1], "%d-%b-%y")

    # def load_caps_data(self):

    # def closest_approach(self):


def cassini_altitude(dt, moon, target='CASSINI'):
    """

    :param dt:
    :param moon:
    :param target:
    :return:
    """
    et = spice.datetime2et(dt)
    frame = 'IAU_' + moon.upper()
    observ = moon.upper()
    corrtn = 'NONE'

    state, ltime = spice.spkpos(target, et, frame, corrtn, observ)
    # TODO do full calculation of altitude
    lon, lat, alt = spice.recpgr(moon.upper(), state, spice.bodvrd(moon.upper(), 'RADII', 3)[1][0], 2.64e-4)

    return alt


def nearest_moon(dt, target='CASSINI'):
    """

    :param dt:
    :param target:
    :return:
    """
    et = spice.datetime2et(dt)
    moondict = {}
    NAIFIDS = range(601,654,1)

    for x in NAIFIDS:
        moon = spice.bodc2s(x)

        frame = 'IAU_' + moon.upper()
        observ = moon.upper()
        corrtn = 'NONE'

        state, ltime = spice.spkpos(target, et, frame, corrtn, observ)
        # TODO do full calculation of altitude
        lon, lat, alt = spice.recpgr(moon.upper(), state, spice.bodvrd(moon.upper(), 'RADII', 3)[1][0], 2.64e-4)
        moondict[moon] = alt
    sorteddict = OrderedDict(sorted(moondict.items(), key=lambda t: t[1]))

    return sorteddict


def moon_alt(start_dt, end_dt, deltat, moon, target='CASSINI'):
    """

    :param end_dt:
    :param start_dt:
    :param deltat:
    :param target:
    :return:
    """
    dtlist = [start_dt]
    altlist = []
    counter = 1
    tempdeltat = deltat
    while start_dt + tempdeltat < end_dt:
        tempdeltat = deltat * counter
        dtlist.append(start_dt + tempdeltat)
        counter += 1

    for tempdt in dtlist:
        altlist.append(cassini_altitude(tempdt, moon))

    return dtlist, altlist


def moon_flybys(moon=""):
    """
    """
    for key, value in flybydict.items():
        if moon == "":
            print(key, value)
        elif value[0] == moon:
            print(key, value)
