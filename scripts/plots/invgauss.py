from scipy.stats import invgauss
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(1, 1)
mu = 0.145462645553
mean, var, skew, kurt = invgauss.stats(scale=2, mu=mu, moments='mvsk')

print (mean, var, skew, kurt)

x = np.linspace(0, 3, 500)
ax.plot(x, invgauss.pdf(x, scale=1,   mu=1), 'blue')
ax.plot(x, invgauss.pdf(x, scale=0.2, mu=1), 'green')
ax.plot(x, invgauss.pdf(x, scale=3,   mu=1), 'red')
ax.plot(x, invgauss.pdf(x, scale=1,   mu=3), 'cyan')
ax.plot(x, invgauss.pdf(x, scale=0.2, mu=3), 'magenta')

plt.show()
