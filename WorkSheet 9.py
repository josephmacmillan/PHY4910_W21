from astropy import coordinates
from astropy.stats import LombScargle
import matplotlib.pyplot as plt
import numpy as np
from astroquery.esa.hubble import ESAHubble

#Imports Hubble data
table =  ESAHubble.query_hst_tap("SELECT * FROM hcv.hcv WHERE matchid=106259089 AND filter LIKE '%ACS_F814W%'")
table.show_in_browser()

#Takes the light curve data form the table
x = table['lightcurve_d']
y = table['lightcurve_cm']
yerr = table['lightcurve_e']

x = x[0:-2]
y= y[0:-2]
yerr = yerr[0:-2]

#Plots the magnitude over period of time
plt.errorbar(x,y,yerr=yerr, marker="o", linestyle="--", linewidth=0.5, color='k')
plt.xlabel("Modified Me(Julian) Date (days)")
plt.ylabel("Modeified Magnitude")
plt.show()

x.unit = None
y.unit = None
yerr.unit = None

#plots the frequency 
ls = LombScargle(x,y,yerr)
freq, power = ls.autopower(minimum_frequency = 1, maximum_frequency = 5)
plt.plot(freq, power)
plt.show()

#Uses the highest frequency to find the period
period = 1/freq[np.argmax(power)]
print(period)

#plots the period
phase = x / period % 1
plt.scatter(phase, y, marker="o", color="k")
plt.scatter(phase+1, y, marker="o", color="k")
plt.show()



