import matplotlib.pyplot as plt
import numpy as np

V = np.load('v.npz')['v']
print(V.shape)
#for i in range(0, 20):
for i in np.random.randint(V.shape[1], size=200):
    #xprint(i)
    plt.plot(V[:, i])
plt.ylim(np.percentile(V, 1), np.percentile(V,90))
print('min,max', V.min(), V.max())
print('median', np.median(V))
print('percentiles', np.percentile(V, [1, 10, 25, 50, 75, 90, 99]))

plt.savefig('sample_v.pdf')


