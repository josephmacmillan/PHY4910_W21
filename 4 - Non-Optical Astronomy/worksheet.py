# CASA commands

# run CASA
LD_PRELOAD=../casa-6.1.1-15-pipeline-2020.1.0.40/lib/libfontconfig.so.1.11.1 ../casa-6.1.1-15-pipeline-2020.1.0.40/bin/casa

# look at channel data 
plotms(vis='Antennae_South.cal.ms',xaxis='channel',yaxis='amp', avgtime='1e8',avgscan=True, showgui = True)

# or show frequency of each channel
plotms(vis='Antennae_South.cal.ms',xaxis='freq',yaxis='amp', avgtime='1e8',avgscan=True, showgui = True)

# remove continuum data ...takes a bit of time
uvcontsub(vis='Antennae_South.cal.ms',fitspw='0:1~30;120~164',fitorder = 1)

# remove old stuff
for ext in ['.image','.model','.image.pbcor','.psf','.residual','.pb','.sumwt','.weight']:
    rmtables('Antennae_South.CO3_2Line.Clean'+ext)

# clean the CO3-2 image data (remember its a cube).  Just do a couple of interactive steps then quit.  
# Takes a long time just to get to the window; went through two full 100%'s.
# clean again, hit green arrow thing: goes through one 100%
tclean(vis='Antennae_South.cal.ms.contsub',
      imagename='Antennae_South.CO3_2Line.Clean',
      spw='0',field='',phasecenter=15,
      specmode='cube',outframe='LSRK',restfreq='345.79599GHz',
      nchan=70,start='1200km/s',width='10km/s',
      gridder='mosaic',deconvolver='hogbom',
      imsize=750,cell='0.13arcsec',pblimit=0.2,
      restoringbeam='common',
      interactive=True,
      weighting='briggs',robust=0.5,
      niter=20000, threshold='5.0mJy',
      savemodel='modelcolumn')

# actually look at the cleaned image?
viewer('Antennae_South.CO3_2Line.Clean.image')
      
# skip the self-calibration ... create the moments, I guess
# then maybe just show the proper results - download them from the website and out them in the worksheet directory?
immoments('Antennae_South.CO3_2Line.Clean.image',
          moments=[0], chans='12~64', 
          outfile='Antennae_South.CO3_2Line.Clean.image.mom0')
          
# view it
viewer('Antennae_South.CO3_2Line.Clean.image.mom0')
          
# save to fits
exportfits(imagename='Antennae_South.CO3_2Line.Clean.image',
     fitsimage='Antennae_South.CO3_2Line.Clean.image.fits')
exportfits(imagename='Antennae_South.CO3_2Line.Clean.image.mom0',
     fitsimage='Antennae_South.CO3_2Line.Clean.image.mom0.fits')
