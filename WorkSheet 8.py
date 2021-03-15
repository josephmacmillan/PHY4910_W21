import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import LogStretch, PercentileInterval
from astropy.stats import sigma_clipped_stats
from photutils import IRAFStarFinder, CircularAperture, aperture_photometry


#Imports fits data
fit = fits.open("hst_12311_02_wfc3_uvis_f275w_drz.fits")
image = fit['SCI'].data

#Looks at a smaller region of image
image = image[3000:4000, 3000:4000]

#filters brightest and darkest 5%
transform = PercentileInterval(99.5)
image_t = transform(image)


mean, median, std =  sigma_clipped_stats(image, sigma=3)

finder = IRAFStarFinder(fwhm = 3, threshold = 5*std)
stars = finder(image - median)
#stars.show_in_browser()


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

#Imports 814w filter fits file
fit = fits.open("hst_12311_02_wfc3_uvis_f814w_drz.fits")
image = fit['SCI'].data

image = image[3000:4000, 3000:4000]

transform = PercentileInterval(99.5)
image_t = transform(image)

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
    
#print(position)

apertures = CircularAperture(position, r=4)
#plt.imshow(image_t, cmap = "gray_r", origin = "lower")
#apertures.plot(color="red", lw = 1.5, alpha = 0.5)

phot=aperture_photometry(image, apertures)
mags_IR=25.0-2.5*np.log10(phot['aperture_sum'])

plt.scatter(mags_UV-mags_IR,mags_UV, marker="o")
plt.gca().invert_yaxis()
plt.show()


