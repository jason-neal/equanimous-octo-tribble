
# Test modify xtick labels
import matplotlib.pyplot as plt
import numpy as np


def convert_labels_to_rv(x, labels):
    x_mid = np.mean(x)
    c = 299792.458   # km/s
    return ["{0:2.1f}".format((label / x_mid) * c) for label in labels]


x = np.arange(10) + 2110
y = (np.arange(10) - 2) * 0.01

ax = plt.subplot(111)

plt.plot(x, y)
ax.set_xlabel("x")
ax.set_ylabel("y")

ax2 = ax.twinx()
ax2.plot(x, y)
ax2.set_ylabel("Transformed (km/s)")

locations = [item for item in ax.get_yticks()]

print(locations)
new_labels = convert_labels_to_rv(x, locations)
ax2.set_yticklabels(new_labels)
plt.show()
