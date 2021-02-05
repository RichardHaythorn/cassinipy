import datetime
import cassinipy.moons as moons
import cassinipy.spice
import spiceypy as spice
import matplotlib.pyplot as plt

startdt = datetime.datetime(2007, 12,3,6)
enddt = datetime.datetime(2007,12,3,8)
dt = datetime.timedelta(minutes=1)

NAIFIDS = range(601, 654, 1)

fig = plt.figure()
ax = fig.add_subplot(111)

for counter, x in enumerate(NAIFIDS):
    print("Moons Checked:", counter)
    moon = spice.bodc2s(x)
    times, alts = moons.moon_alt(startdt, enddt, dt, moon)
    if min(alts) < 1e4 and moon != "TITAN":
        ax.plot(times, alts, label=moon)

ax.set_yscale("log")
ax.grid()
ax.legend()
plt.show()
