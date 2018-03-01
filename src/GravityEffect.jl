module GravityEffect

import FileTools
import DataFrames
import ResampleAndFit


include("eopeffects.jl");
include("atmacseffect.jl");
include("eosteffects.jl");
include("inpolygon.jl");
include("tesseroid.jl");
include("simplebodies.jl");

export prismEffect, cylinderEffect, bouguerEffect, pointEffect, distance,
	   prism2point, correctCurvature, curvatureEffect, tesseroid,
	   eopEffect,
	   atmacsEffect,
	   eostEffect,
	   inpolygon

# Earth's radius and Gravitational constant
const R_const = 6371000.; # m
const G_const = 6.674215*10^(-11.); # Nm^2/kg^2
const W_const = 72921151.467064/10^12; # angular velocity

end # module
