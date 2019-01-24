GravtityEffect
================================================================
[![Build Status](https://travis-ci.org/emenems/GravityEffect.jl.svg?branch=master)](https://travis-ci.org/emenems/GravityEffect.jl)
[![codecov](https://codecov.io/gh/emenems/GravityEffect.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/emenems/GravityEffect.jl)
[![Coverage Status](https://coveralls.io/repos/github/emenems/GravityEffect.jl/badge.svg?branch=master)](https://coveralls.io/github/emenems/GravityEffect.jl?branch=master)


This Julia package contains functions for computation of gravity effect of
following bodies:

* Bouguer plate: `bouguerEffect`
* Cylinder: `cylinderEffect`
* Point mass (or sphere): `pointEffect`
* Prism (in Cartesian coordinate system) `prismEffect`
* Tesseroid in geodetic (longitude, latitude, height) coordinates `tesseroid`  

To compute gravity effect time series use:
* `layerResponse`: gravity response to (soil) layers copying the terrain (digital
elevation model) in arbitrary depths. The output of this function can be used in:
* `sm2effect`: convert soil moisture (m<sup>3</sup>/m<sup>3</sup>) to gravity effect time series (nm/s<sup>2</sup>)  

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
* aggregate unit `layerResponse` to required depth range: `aggregate_layers`

In addition, following constants are used (not global):
* `G`: gravity constant
* `R`: radius of Earth replacement sphere

## Usage
* Check the function help for instructions and example usage, e.g., `?bouguerEffect`

> Check the `REQUIRE` file for required packages  
> [FileTools.jl](https://github.com/emenems/FileTools.jl) and [ResampleAndFit.jl](https://github.com/emenems/ResampleAndFit.jl) packages can downloaded from GitHub
