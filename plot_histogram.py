import matplotlib.pyplot as plt
import numpy as np

f = open("out.txt", "r")
bin_size = float(f.readline())
print(bin_size)
log_freqs = []
for s in f.readline().split():
    log_freqs.append(float(s))

condVals = bin_size*np.array(range(len(log_freqs)))
# #*bin_size
# ax, fig = plt.subplots(111)
plt.plot(condVals,np.exp(log_freqs))
plt.xlabel(f'Conductance (bin size {bin_size})')
plt.ylabel('log(# partitions)')
plt.yscale('log')
plt.xlim([0,1])
plt.ylim([0.5,10**6])
plt.savefig('log_freqs.png')