import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import LogStretch, PercentileInterval
from astropy.stats import sigma_clipped_stats
from photutils import IRAFStarFinder, CircularAperture, aperture_photometry


"""
Importing fits data
"""
fit = fits.open("hst_12311_02_wfc3_uvis_f275w_drz.fits")
image = fit['SCI'].data

"""
slicing the array to deal with a smaller image
"""
image = image[3000:4000, 3000:4000]


"""
normalization of stretch of the image
"""
transform = PercentileInterval(99.5)
image_t = transform(image)

"""
using a method called sigma clipping to get statistics about the background of an image
"""
mean, median, std =  sigma_clipped_stats(image, sigma=3)

finder = IRAFStarFinder(fwhm = 3, threshold = 5*std)
stars = finder(image - median)
#stars.show_in_browser() #this shows a whole array


positions = np.array((stars["xcentroid"],stars["ycentroid"]))

position = []
for i in range(len(stars)):
    x = stars["xcentroid"][i]
    y = stars["ycentroid"][i]
    position.append((x,y))
    
#print(position)

#Circle the found stars
apertures = CircularAperture(position, r=4)
#plt.imshow(image_t, cmap = "gray_r", origin = "lower")
#apertures.plot(color="red", lw = 1.5, alpha = 0.5)

phot=aperture_photometry(image, apertures)
mags_UV=25.0-2.5*np.log10(phot['aperture_sum'])

#################################################

"""
Importing fits data
"""
fit = fits.open("hst_12311_02_wfc3_uvis_f814w_drz.fits")
image = fit['SCI'].data

"""
splicing, same as above
"""
image = image[3000:4000, 3000:4000]

"""
normalization of stretch of the image
"""

transform = PercentileInterval(99.5)
image_t = transform(image)

"""
this bit to show the sigma clipped data
"""

#mean, median, std =  sigma_clipped_stats(image, sigma=3)

#finder = IRAFStarFinder(fwhm = 3, threshold = 5*std)
#stars = finder(image - median)
#stars.show_in_browser()

positions = np.array((stars["xcentroid"],stars["ycentroid"]))

position = []
for i in range(len(stars)):
    x = stars["xcentroid"][i]
    y = stars["ycentroid"][i]
    position.append((x,y))
    
#print(position) #just to check

apertures = CircularAperture(position, r=4)
#plt.imshow(image_t, cmap = "gray_r", origin = "lower")
#apertures.plot(color="red", lw = 1.5, alpha = 0.5)

"""
now we perform photometry to find the flux of each star from the array, essentially getting the sum of all pixel values within the aperture 
"""
phot=aperture_photometry(image, apertures)
mags_IR=25.0-2.5*np.log10(phot['aperture_sum'])

plt.scatter(mags_UV-mags_IR,mags_UV, marker="o")
plt.gca().invert_yaxis()
plt.show()