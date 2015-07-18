import matplotlib.pyplot as plt
import numpy as np

V = np.load('sample_v.npy')
for i in range(0, 100):
    plt.plot(V[i, :])
plt.ylim(0, 100)
print(V)
plt.savefig('sample_v.pdf')


