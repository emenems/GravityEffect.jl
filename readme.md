GravtityEffect: Julia package for gravity effects
================================================================
This repository contains functions for computation of gravity effect of following bodies:

* Bouguer plate: `bouguerEffect`
* Cylinder: `cylinderEffect`
* Point mass (or sphere): `pointEffect`
* Prism (in Cartesian coordinate system) `prismEffect`

> Check the function help for instructions and example usage

Following auxiliary functions are provided:
* distance between points (2D or 3D): `distance`
* convert prism to point mass (convert density and prism resolution to mass) `prism2point`
* compute Earth's curvature effect on height (Cartesian coordinates are used): `curvatureEffect`
* correct height for Earth's curvature: `correctCurvature`

In addition, following constants are used (not global):
* `G`: gravity constant
* `R`: radius of Earth replacement sphere
