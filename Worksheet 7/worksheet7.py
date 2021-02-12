import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.visualization import LogStretch, PercentileInterval, AsinhStretch

# Let's open up the FITS file:
f = fits.open("hst_12543_01_wfc3_uvis_f658n_sci.fits")

# Now, what exactly is in this FITS file?  Check it out:
print(f.info())

# Okay, looks like two HDUs, one called PRIMARY (with no image data!) and one called SCI (with image data).  Load the SCI:
image  = f['SCI'].data
header = f['SCI'].header

# The header contains lots of cool info:
print(header)

# We can look at the image with MatPlotLib.  But it's better to apply a transformation to the
# image first; take a look at the AstroPy docs to see what transforms we can use.

transform = AsinhStretch(0.9) + PercentileInterval(99.5)

# Now show the image.  We'll read the coordinate system info from the header.
# And I'll choose a reverse grayscale colormap and set the origin to be the bottom left corner.
ax = plt.subplot(projection = WCS(header))
ax.imshow(transform(image), cmap='gray_r', origin='lower')
plt.show()


