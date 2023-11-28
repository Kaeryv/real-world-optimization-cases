# cec_2011_real_world_cases

# Be careful when using: Work in progress

Only functions 1 to 7 are available.

## Installation

```bash
pip install 'cec_2011_rw @ git+https://github.com/Kaeryv/cec_2011_real_world_cases'
```

## Usage

```python
from cec_2011_rw.real_world_cases import function_selector, get_boundaries
import numpy as np
import matplotlib.pyplot as plt


# Test varying the first parameter of each problem.
fig, axs = plt.subplots(2, 4)
axs = axs.flatten()
for i in range(1, 8):
    bounds = get_boundaries(i)
    d = bounds.shape[1]
    fs = list()
    x0 = np.linspace(bounds[0, 0], bounds[1, 0], 256)
    for x in x0:
        params = bounds[0] + np.random.rand((d)) * (bounds[1] - bounds[0])
        params[0] = x
        f = function_selector(params, num=i)
        fs.append(f)
    axs[i-1].plot(x0, fs)
    axs[i-1].set_title(f"rw_{i}")
plt.tight_layout()
fig.savefig("test_x0.png")
```
