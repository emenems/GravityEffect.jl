GravtityEffect: Julia package for gravity effects
================================================================
This repository contains functions for computation of gravity effect of following bodies:

* Bouguer plate: `bouguerEffect`
* Cylinder: `cylinderEffect`
* Point mass (or sphere): `pointEffect`
* Prism (in Cartesian coordinate system) `prismEffect`
* Tesseroid in geodetic (longitude, latitude, height) coordinates `tesseroid`  

A complete gravity response to (soil) layers copying the terrain (digital elevation model) in arbitrary depths can be computed using `layerResponse` function.

The repository also contains files for computation/downloading of:

* Polar motion and Length of Day effects: `eopEffects`
* [Atmacs](http://atmacs.bkg.bund.de/docs/data.php) atmospheric effect on surface gravity: `atmacsEffect`
* [EOST Loading](http://loading.u-strasbg.fr/surface_gravity.php) effects (atmosphere, hydrology, non-tidal ocean loading and pol/LOD effects)

> Check the function help for instructions and example usage

Following auxiliary functions are provided:
* distance between points (2D or 3D): `distance`
* convert prism to point mass (convert density and prism resolution to mass) `prism2point`
* compute Earth's curvature effect on height (Cartesian coordinates are used): `curvatureEffect`
* correct height for Earth's curvature: `correctCurvature`
* find points inside a polygon: `inpolygon`

In addition, following constants are used (not global):
* `G`: gravity constant
* `R`: radius of Earth replacement sphere

> Check the `REQUIRE` file for required packages
