# ZonoOpt
This C++ library provides classes and tailored optimization routines for zonotopes, constrained zonotopes, and hybrid zonotopes. 

![image](images/mhe-traj.svg)
![image](images/reachability-traj.svg)

## Zonotopes, Constrained Zonotopes, and Hybrid Zonotopes
Zonotopes and generalizations of zonotopes are set representations widely used for reachability analysis. 
This library focuses specifically on zonotopes, constrained zonotopes, and hybrid zonotopes. 
A set $\mathcal{Z}$ is a zonotope if there exists a generator matrix $G$ and center vector $c$ such that
$$ \mathcal{Z} = \left\{ G \xi + c \middle| \xi \in [-1, 1]^{n_G} \right\} \;. $$

TODO: constrained zonotopes, hybrid zonotopes, images illustrating each

## ZonoOpt Features

All classes and methods are implemented using sparse linear algebra via the Eigen library.

ZonoOpt uses polymorphism to provide a common interface for zonotopes, constrained zonotopes, and hybrid zonotopes.

Generators may optionally have range $[0,1]$ instead of $[-1,1]$.

ZonoOpt has no external dependencies beyond Eigen, making it easy to integrate into robotics projects using C++ or Python.

## Building and Installing
Python bindings can be installed from PyPI with `pip install zonoopt`. To build the bindings from source, use `pip install .`. Note that a C++ compiler is required to build from source.

This library can be used in CMake projects either via add_subdirectory or by installing the library. Including the library via add_subdirectory is recommended when possible as it permits more aggressive CPU optimizations (i.e., -march=native). To link your application to ZonoOpt, use `target_link_libraries(your_application PRIVATE ZonoOpt)`. When building the library for installation, you must set the option ZONOOPT_INSTALL to ON, i.e., `cmake -DZONOOPT_INSTALL=ON -S . -B build`. 

More information about ZonoOpt can be found in the following publication. Please cite this if you publish work based on ZonoOpt: 
**Robbins, J.A., Siefert, J.A., and Pangborn, H.C., "Sparsity-Promoting Reachability Analysis and Optimization of Constrained Zonotopes," 2025. [https://arxiv.org/abs/2504.03885](https://doi.org/10.48550/arXiv.2504.03885).**

## Examples
Consider the case that we wish to compute the robust forward reachable set of a discrete time double integrator system and verify that it does not intersect an unsafe set.
We may do this in Python as follows:
```python
import zonoopt as zono
import numpy as np
import matplotlib.pyplot as plt

# System dynamics
dt = 0.1
A = np.array([[1., dt],
              [0., 1.]])
B = np.array([[0.5*dt**2],
              [dt]])

# Initial set: box [-1.0, 1.0] x [-0.1, 0.1]
X0 = zono.interval_2_zono(zono.Box([-1., -0.1], [1., 0.1]))

# Input set: box [-0.2, 0.2]
U = zono.interval_2_zono(zono.Box([-0.2], [0.2]))

# Disturbance set: affine map of octagon with radius 0.05
W = zono.make_regular_zono_2D(radius=1., n_sides=8)
W = zono.affine_map(W, np.diag([0.01, 0.05]))

# Compute reachable set over 10 time steps
X = X0
for k in range(10):
    X = zono.affine_map(X, A)
    X = zono.minkowski_sum(X, zono.affine_map(U, B))
    X = zono.minkowski_sum(X, W)

# Unsafe set
O = zono.vrep_2_conzono(np.array([[1.3, 0.],
                                  [1.6, 0.8],
                                  [2.0, -0.4],
                                  [2.3, 0.6]]))

# Check for intersection with unsafe set
print(f'10-step reachable set intersects unsafe set: {not zono.intersection(X, O).is_empty()}')

# Plot the final reachable set
zono.plot(X, color='b', alpha=0.2)
zono.plot(O, color='r', alpha=0.2)
plt.show()
```

TODO: Add C++ example

Extensive Python examples are located in examples, and a C++ example is located in test.

## Documentation
Auto-generated API documentation is available below.

[C++ API](https://psu-PAC-Lab.github.io/ZonoOpt/C++/html/index.html)

[Python API](https://psu-PAC-Lab.github.io/ZonoOpt/python/build/html/index.html)
