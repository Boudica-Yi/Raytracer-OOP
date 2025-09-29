# Optical Ray Tracer

This repository implements a Python-based geometrical ray tracer to simulate the propagation of light rays through different optical elements and to investigate key aberrations such as spherical aberration and coma. The project includes methods for propagation of single rays, generation of collimated ray bundles, focal point calculation, RMS spot radius and diffraction limit analysis, and lens curvature optimisation for improved imaging quality.

---

## Project Overview

The raytracer module computes the trajectories of rays through optical systems composed of spherical surfaces and plano-convex or biconvex lenses. Key implemented features include:

- Ray propagation through refracting surfaces using Snell’s Law.  
- Uniform collimated ray bundle generation with circular sampling.  
- Focal length calculation using both Lens Maker’s formula and graphical methods.  
- RMS spot radius and diffraction limit calculations for evaluating image quality.  
- Simulation of monochromatic aberrations including spherical and coma aberrations.  
- Lens curvature optimisation by minimising RMS spot size across curvature pairs.

The tracer is written entirely in Python and uses object-oriented programming for lens inheritance and special cases (e.g. plano-convex reversed orientation).

---

## Project Structure

```text
.
├── raytracer.py              # Core raytracer module (all classes and functions)
├── Tesing2_T11_T14.py        # Tests for Tasks 11–14
├── Tesing3_Task15.py         # Tests for Task 15
├── Testing4.py               # Lens optimisation tests (long runtime)
├── Raytracer Plots           # Graphical illustration of all plots from tests
└── README.md                 # This file

```

---

## Methods Implemented

### 1. Single Ray Propagation
Rays are represented by positions and directions. At each refracting surface, the intersection point and new direction are calculated using Snell’s law:

$$
\vec{k}_r = \frac{n_2}{n_1}\vec{k}_i + \left(\frac{n_2}{n_1}\cos\theta_1 - \cos\theta_2\right)\vec{n}
$$

where \( n_1, n_2 \) are refractive indices, \(\vec{n}\) is the surface normal, and \(\theta_1, \theta_2\) are incident and refraction angles.

---

### 2. Ray Bundles
Collimated ray bundles are generated using layered circular sampling, ensuring uniform spatial distribution across the bundle cross-section. Each bundle ray inherits the central axis direction for propagation.

---

### 3. Focal Length Estimation
Two methods are implemented:

- Analytical: Lens Maker’s Formula, used for initial estimates.  
- Graphical: Propagating a near-paraxial ray and finding its axis intercept. This method is used throughout due to spherical aberration effects.

---

### 4. Aberration Analysis
The module simulates and visualises:

- Spherical aberration: rays far from axis focus at different points.  
- Coma aberration: off-axis bundles produce asymmetric foci.

For example, rays at 5 mm bundle diameter and 700 nm wavelength show significant paraxial deviation:contentReference[oaicite:1]{index=1}.

---

### 5. RMS Spot Radius and Diffraction Limit
The RMS spot radius of the bundle at the focal plane is compared against the diffraction limit:

$$
L_\text{diff} = \frac{\lambda f}{D}
$$

where \( \lambda \) is wavelength, \( f \) focal length, and \( D \) bundle diameter.

Example values:contentReference[oaicite:2]{index=2}:

| Wavelength | r<sub>RMS</sub> (mm) | L<sub>diff</sub> (mm) |
|-----------|----------------------|-----------------------|
| Red (700 nm) | 0.003 | 0.014 |
| Blue (445 nm) | 0.002 | 0.009 |

---

### 6. Lens Optimisation
For biconvex lenses, the code minimises r<sub>RMS</sub> over curvature pairs via:

- Method 1: Numerical minimisation  
- Method 2: Parameter sweep

Both methods return similar optimal curvatures (e.g. r<sub>min</sub> ≈ 0.0025 mm at focal plane z ≈ 120 mm):contentReference[oaicite:3]{index=3}.

