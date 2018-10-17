import numpy as np
import matplotlib.pyplot as plt
x = np.linspace(-np.pi, np.pi, 101)
y1 = np.cos(x)  # Numpy math
y2 = np.sin(x)  # operate on whole
plt.plot(x, y1, color="r", linestyle="-", linewidth=2)
plt.plot(x, y2, "g--", lw=1)
# shorter description: thinner green dashed
plt.show()