from pycbc import catalog
import matplotlib.pyplot as plt

# this worksheet is for going over the basics of signal analysis

m = catalog.Merger('GW150914')

h = m.strain('H1')
plt.plot(h.sample_times, h)
plt.show()

# zoom into the one second around the merger time
hz = h.time_slice(m.time - 0.5, m.time + 0.5)
plt.plot(hz.sample_times, hz)
plt.show()

# look at the power spectral density
psd = h.psd(4)
plt.loglog(psd.sample_frequencies, psd)
plt.show()

# whiten the data ...
hw = h.whiten(4, 4)  # 4 seconds of each sample used in PSD estimate, 4 second duration of filter
psdw = hw.psd(4)
#plt.loglog(psd.sample_frequencies, psd)
plt.loglog(psdw.sample_frequencies, psdw)
plt.show()

# what does the strain look like now?
#plt.plot(hw.sample_times, hw, linewidth=0.5)
hwz = hw.time_slice(m.time - 0.5, m.time + 0.5)
plt.plot(hwz.sample_times, hwz)
plt.show()

# low and high pass filters
hf = hw.highpass_fir(30, 512).lowpass_fir(250, 512)
hfz = hf.time_slice(m.time - 0.5, m.time + 0.5)

plt.plot(hfz.sample_times, hfz)
plt.show()

# check by doing the same thing for L1:
hl = m.strain('L1')
h = h.whiten(4,4).highpass_fir(30, 512).lowpass_fir(250, 512).time_slice(m.time - 0.5, m.time + 0.5)
hl.roll(int(hl.sample_rate * .007))
hl *= -1
plt.plot(hfz.sample_times, hfz)
plt.plot(hl.sample_times, hl)
plt.show()


# finally, Q-transform plot
hfz =  hw.time_slice(m.time - 0.5, m.time + 0.5)
times, freqs, power = hfz.qtransform(.001, logfsteps=100,
                                                qrange=(8, 8),
                                                frange=(20, 512),
                                         )
plt.pcolormesh(times, freqs, power**0.5)
plt.yscale('log')
plt.show()
