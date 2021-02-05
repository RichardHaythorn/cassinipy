import numpy as np
import spiceypy as spice


def rotate_CAPS_SCframe(CAPS_actuation, instrument, anodes=False):
    """
    Returns the CAPS pointing vector after actuation in the SC frame
    """
    naifiddict = {'ims': -82820, 'els': -82821, 'ibs1': -82822, 'ibs2': -82823, 'ibs3': -82824}
    naifid = naifiddict[instrument]

    if not anodes:
        rotationmatrix = spice.spiceypy.axisar(np.array([0, 0, -1]), CAPS_actuation * spice.rpd())
        temp = spice.mxv(rotationmatrix, spice.spiceypy.getfov(naifid, 20)[2])

    if (instrument == 'ims' or instrument == 'els') and anodes == True:
        temp = []

        rotationmatrix_act = spice.spiceypy.axisar(np.array([0, 0, -1]), CAPS_actuation * spice.rpd())

        for anodenumber, x in enumerate(np.arange(70, -90, -20)):
            rotationmatrix_anode = spice.spiceypy.axisar(np.array([-1, 0, 0]), x * spice.rpd())
            postanode_rotation = spice.mxv(rotationmatrix_anode, spice.spiceypy.getfov(naifid, 20)[2])

            temp.append(spice.mxv(rotationmatrix_act, postanode_rotation))

    return temp


def cassini_ramdirection_SCframe(tempdatetime, target='CASSINI', frame='J2000', observ='titan', corrtn='NONE',
                                 output=False):
    et = spice.datetime2et(tempdatetime)

    state, ltime = spice.spkezr(target, et, frame, corrtn, observ)
    ramdir = spice.vhat(state[3:6])

    # Gets Attitude
    sclkdp = spice.sce2c(-82, et)  # converts an et to a continuous encoded sc clock (ticks)
    ckgp_output = spice.ckgp(-82000, sclkdp, 0, frame)
    cmat = ckgp_output[0]

    ram_unit = spice.mxv(cmat, ramdir)

    if output:
        print('ET        = {:20.6f}'.format(et))
        print('VX        = {:20.6f}'.format(ram_unit[0]))
        print('VY        = {:20.6f}'.format(ram_unit[1]))
        print('VZ        = {:20.6f}'.format(ram_unit[2]))

    return ram_unit


def caps_ramdirection_azielv(tempdatetime, observ='titan'):
    """
    Returns the azimuth and elevation of the ramdirection
    """
    ram_unit = cassini_ramdirection_SCframe(tempdatetime, observ=observ)

    sc2CAPS = np.array(([0, -1, 0],
                        [-1, 0, 0],
                        [0, 0, -1]))
    ram_unit_CAPS = spice.mxv(sc2CAPS, ram_unit)
    # print(ram_unit_CAPS)
    ram_unit_azielv = np.array(spice.reclat(ram_unit_CAPS)[1:]) * spice.dpr()

    return ram_unit_azielv[0], ram_unit_azielv[1]

