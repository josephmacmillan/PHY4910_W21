import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import LogStretch, PercentileInterval, AsinhStretch, make_lupton_rgb

files = ["hst_11577_20_wfc3_uvis_f555w_sci.fits", "hst_11577_20_wfc3_uvis_f625w_sci.fits", "hst_11577_20_wfc3_uvis_f814w_sci.fits"]

transform = AsinhStretch(0.7) + PercentileInterval(99.5)

images = []
master_header = None
for i in range(len(files)):
    print("Opening FITS file", files[i])
    with fits.open(files[i]) as hdul:
        hdu = hdul[1]
        if i == 0:
            print("Setting master image and header", i)
            master_header = hdu.header
        image = hdu.data
        images.append(transform(image))

rgb = make_lupton_rgb(images[0], images[1], images[2], Q = 10, minimum = 0.0)
fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection=WCS(master_header))
ax.imshow(transform(rgb), origin='lower', interpolation='bilinear')

hdul = fits.open("Antennae_North.CO3_2Line.Clean.pcal1.image.mom0.fits")
radio = hdul[0].data[0][0]
radio_header = hdul[0].header
radio_wcs = WCS(radio_header)
radio_wcs = radio_wcs.sub(2)
ax.contour(radio, origin='lower', colors="black", levels=[-9.19579, 5.68717, 20.5701, 35.4531, 50.336, 65.219, 80.1019, 94.9849], transform=ax.get_transform(radio_wcs))

hdul = fits.open("Antennae_South.CO3_2Line.Clean.pcal1.image.mom0.fits")
radio = hdul[0].data[0][0]
radio_header = hdul[0].header
radio_wcs = WCS(radio_header)
radio_wcs = radio_wcs.sub(2)
ax.contour(radio, origin='lower', colors="black", levels=[-9.19579, 5.68717, 20.5701, 35.4531, 50.336, 65.219, 80.1019, 94.9849], transform=ax.get_transform(radio_wcs))
#ax.imshow(radio, origin='lower', interpolation='bilinear', transform=ax.get_transform(radio_wcs))

plt.show()


