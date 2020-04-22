import spiceypy as spice

def cassini_phase(utc,target='CASSINI',frame='IAU_TITAN',observ='TITAN',corrtn='NONE',output=False):
     
    et = spice.str2et(utc)    
    
    state, ltime = spice.spkezr(target, et, frame,corrtn, observ)
    
    if output == True:
        print('ET        = {:20.6f}'.format(et))
        print(' X        = {:20.6f}'.format(state[0]))
        print(' Y        = {:20.6f}'.format(state[1]))
        print(' Z        = {:20.6f}'.format(state[2]))
        print('VX        = {:20.6f}'.format(state[3]))
        print('VY        = {:20.6f}'.format(state[4]))
        print('VZ        = {:20.6f}'.format(state[5]))
        
    return state

def cassini_surfintercerpt(utc,output=False):
    
    target = 'TITAN'
    fixref  = 'IAU_TITAN'         
    dref  = 'IAU_TITAN'        
    method = 'ELLIPSOID'
    abcorr = 'NONE'
    obsrvr = 'CASSINI'
    state = cassini_phase(utc)
    dvec = spice.vhat(-state[:3])
   
    et = spice.str2et(utc)

    [point, trgepc, srfvec] = spice.sincpt(method, target, et, fixref, abcorr, obsrvr, dref, dvec)   
    temp = spice.reclat(point)
    
    radius = temp[0]
    colat = 90+spice.convrt(temp[1], "RADIANS", "DEGREES")
    lon = np.mod(spice.convrt(temp[2], "RADIANS", "DEGREES"),360)
    
    if output == True:
        print("Distance to Intercept Point", spice.vnorm(srfvec))
        print("Radius", radius)
        print("Colatitude", colat )
        print("Longtitude", lon)

    return point, [radius,colat,lon], spice.vnorm(srfvec)

def cassini_altlatlon(utc,output=False):
    
    target = 'TITAN'
    fixref  = 'IAU_TITAN'         
    dref  = 'IAU_TITAN'        
    method = 'ELLIPSOID'
    abcorr = 'NONE'
    obsrvr = 'CASSINI'
    state = cassini_phase(utc)
    
    lon, lat, alt = spice.recpgr(target,state[:3],spice.bodvrd('titan','RADII',3)[1][0],2.64e-4)
    
    if output == True:
        print("ALtitude", alt)
        print("Latitude", lat*spice.dpr())
        print("Longtitude", lon*spice.dpr())

    return alt, lat*spice.dpr(), lon*spice.dpr()
    