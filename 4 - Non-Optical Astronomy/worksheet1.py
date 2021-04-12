from pycbc import catalog
import matplotlib.pyplot as plt

# print a list of all the available catalogs ...

for c in catalog.catalog.list_catalogs():
	print(c)

for merger_name in catalog.Catalog(source='gwtc-1'):
    print(merger_name)
    
for merger_name in catalog.Catalog(source='gwtc-2'):
    print(merger_name)
    
# get some overall data ...
c = catalog.Catalog(source='gwtc-2')
mchirp = c.median1d('mchirp')
dist = c.median1d('distance')
spin = c.median1d('chi_eff')
plt.scatter(dist, mchirp)
plt.ylabel("Chirp Mass")
plt.xlabel("Distance")
plt.show()

# get some specific data - let's do that giant BH GW190521-v3
merger = catalog.Merger("GW190521-v3", source='gwtc-2')
m1 = merger.median1d('mass1')
m2 = merger.median1d('mass2')
print(m1, m2, merger.time)

# get the strain from the Hanford observatory
ts = merger.strain('H1').time_slice(merger.time - 10, merger.time + 10)
print(ts)

print('Duration={:.0f}s Interval={:.1e}s '.format(ts.duration, ts.delta_t))

plt.plot(ts.sample_times, ts)
plt.show()

# now let's look at that neutron star MMM:
merger = catalog.Merger("GW170817-v3", source='gwtc-1')
m1 = merger.median1d('mass1')
m2 = merger.median1d('mass2')
print(m1, m2, merger.time)

# get the strain from the Hanford observatory
ts = merger.strain('H1').time_slice(merger.time - 10, merger.time + 10)
print(ts)

print('Duration={:.0f}s Interval={:.1e}s '.format(ts.duration, ts.delta_t))

plt.plot(ts.sample_times, ts)
plt.show()


# next time - signal analysis

