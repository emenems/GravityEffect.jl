module GravityEffect

import FileTools
import DataFrames
import ResampleAndFit

include("eopeffects.jl");
include("atmacseffect.jl");

export prismEffect, cylinderEffect, bouguerEffect, pointEffect, distance,
	   prism2point, correctCurvature, curvatureEffect, tesseroid,
	   eopEffect,
	   atmacsEffect

# Earth's radius and Gravitational constant
const R_const = 6371000.; # m
const G_const = 6.674215*10^(-11.); # Nm^2/kg^2
const W_const = 72921151.467064/10^12; # angular velocity

include("tesseroid.jl");

"""
	cylinderEffect(depth,radius,thick,density)

Compute gravity effect of a cylinder in m/s^2

**Input**
* depth: depth of the cylinder below the sensor = distance between sensor and the upper boundary in meters
* radius: radius of the cylinder in meters
* thick: cylinder thickness in meters
* density: differential density of the prism in kg/m^3
* dg: gravity effect in m/s^2

**Output**
* gravity effect in m/s^2

**Example**

```
# Compute effect for 1 m depth, radius 100 m, thickness of 1 m and density 1000 kg/m^3
dg = cylinderEffect(1.,100.,1.,1000.)*1e+9;

```
"""
function cylinderEffect{T<:Float64}(depth::T,radius::T,thick::T,density::T)
    return 2*pi*G_const*density*(thick + sqrt(depth^2 + radius^2) - sqrt((depth + thick)^2 + radius^2));
end


"""
	bouguerEffect(thick,density)

Compute Bouguer plate effect

**Input**
* thick: plate thickness in meters
* density: plate density in kg/m^3

**Output**
* dg: gravity effect in m/s^2

**Example**

```
dg = bouguerEffect(1.,1000.)*1e+9;
```
"""
function bouguerEffect{T<:Float64}(thick::T,density::T)
    return 2*pi*G_const*density*thick;
end

"""
	prismEffect(point,prism,resol,density)

Compute gravity effect of a prism in Cartesian coordinate system

**Input**
* sensor: (x,y,z) coordinates of the sensor in meters
* prism: (x,y,z): x and y coordinate of the center of the prism in meters, z coordinate points to upper bound not center!! All in meters.
* resol: (dx,dy,dz) size of the prism (resolution) in x,y,z directions in meters. dz gives the depth of the prism (all in meters)
* density: differential density of the prism in kg/m^3 (defalt value=1.0)

**Output**
* dg: gravity effect in m/s^2

**Example**
```
pe = prismEffect((0.,0.,1.), (0.,0.,0.), (10.,10.,2.), 1000.)*1e+9

```
"""
function prismEffect(sensor::Tuple{Float64,Float64,Float64},
					 prism::Tuple{Float64,Float64,Float64},
					 resol::Tuple{Float64,Float64,Float64},
					 density::Float64=1.0)
    # Shift coordinates system towards the point of computation
    zi = -prism[3] + sensor[3];
    xi = prism[1] - sensor[1];
    yi = prism[2] - sensor[2];
    # Call function for computing gravity effect according to Sorokin formula.
    # Input grids are shifter (+/- grid/2) to point to the edges of the prism
    return density*sorokin(yi-resol[2]/2.,xi-resol[1]/2.,zi,yi+resol[2]/2.,xi+resol[1]/2.,zi+resol[3]);
end

"""
Sorokin's formula for gravity effect of a prism (aux. function)
"""
function sorokin(x1::Float64,y1::Float64,z1::Float64,x2::Float64,y2::Float64,z2::Float64)
    g_log = x2*(log(prismLog(y2,x2,z2)/prismLog(y2,x2,z1)) - log(prismLog(y1,x2,z2)/prismLog(y1,x2,z1))) -
            x1*(log(prismLog(y2,x1,z2)/prismLog(y2,x1,z1)) - log(prismLog(y1,x1,z2)/prismLog(y1,x1,z1))) +
            y2*(log(prismLog(x2,y2,z2)/prismLog(x2,y2,z1)) - log(prismLog(x1,y2,z2)/prismLog(x1,y2,z1))) -
            y1*(log(prismLog(x2,y1,z2)/prismLog(x2,y1,z1)) - log(prismLog(x1,y1,z2)/prismLog(x1,y1,z1)));
    g_tan = z2*(prismAtan(z2,x2,y2)-prismAtan(z2,x1,y2)-prismAtan(z2,x2,y1)+prismAtan(z2,x1,y1)) -
            z1*(prismAtan(z1,x1,y1)-prismAtan(z1,x1,y2)-prismAtan(z1,x2,y1)+prismAtan(z1,x2,y2));
    return -G_const*(g_log + g_tan);
end

"""
	pointEffect(sensor,point,mass)

Compute gravity effect of a point mass

**Input**
* sensor: (x,y,z) coordinate of the sensor in meters
* point: (x,y,z) coordinate of the point mass in meters
* mass: total mass in kg

**Output**
* dg: gravity effect in m/s^2

**Example**:
```
dg = pointEffect((0.,0.,0.), (0.,1.,0.), 1000.)*1e+9;
```
"""
function pointEffect(sensor::Tuple{Float64,Float64,Float64},
					 point::Tuple{Float64,Float64,Float64},
					 mass::Float64=1.0)
    dist = distance(sensor,point);
    dg0 = (mass*G_const)/(dist^2); # direct effect (not in vertical direction)
    cos_alpha = (sensor[3]-point[3])/(dist); # angle between sensor and mass
    return dg0.*cos_alpha;
end

"""
	prism2point(resol,density)

Convert prism to point mass.
**Warning**: to accuratelly compute the gravity effect using point approximation,
correct the height for vertical resolutin (`point[3] = prism[3]-resol[3]/2`)!

**Input**
* resol: (dx,dy,dz) size of the prism (resolution) in x,y,z directions in meters
* density: differential density of the prism in kg/m^3 (defalt value=1.0)

**Output**
* mass in kg

**Example**
```
m = prism2point((1.,1.,1.),10.);
```
"""
function prism2point(resol::Tuple{Float64,Float64,Float64},
					 density::Float64=1.0)
    return mass = resol[1]*resol[2]*resol[3].*density;
end

"""
	distance(point1,point2)

Compute planar or spatial distance between two points

**Input**
* point1: (x,y,z) coordinates of the first point in meters. z coordinate is optional (=> spatial).
* point2: (x,y,z) coordinates of the second point in meters. z coordinate is optional (=> spatial).

**Output**
* planar (2D) or spatial (3D) distnace in meters

**Example**
```
d2 = distance((0.,0.),(0.,1.))
d3 = distance((0.,0.,0.),(0.,0.,1.))
```
"""
function distance(point1::Tuple{Float64,Float64},point2::Tuple{Float64,Float64})
    return sqrt((point1[1]-point2[1])^2 + (point1[2]-point2[2])^2);
end
function distance(point1::Tuple{Float64,Float64,Float64},point2::Tuple{Float64,Float64,Float64})
    return sqrt((point1[1]-point2[1])^2 + (point1[2]-point2[2])^2 + (point1[3]-point2[3])^2);
end

"""
	curvatureEffect(dist)

Effect of Earth's curvature on height as a function of distance

**Input**
* dist: distance to the point where height should be corrected (in meters)

**Output**
* curvatrure effect (change of height) in meters

**Example**
```
ce = curvatureEffect(10000.0)
```
"""
function curvatureEffect(dist::Float64)
    return ((sqrt(dist^2 + R_const^2) - R_const).*R_const)/(sqrt(dist^2 + R_const^2));
end

"""
	correctCurvature(sensor,grid)

Correct height for Earth's curvature

**Input**
* sensor: (x,y,z) coordinate of the sensor in meters (z is optional and will not be used)
* grid: (x,y,z) coordinate of the point mass or prism center in meters

**Output**
* corrected height (grid z coordinate) in meters

**Example**
```
hcorr = correctHeight((0.,0.),(0.,10000.,0.))
hcorr = correctHeight((0.,0.,0.),(0.,10000.,0.))
```
"""
function correctCurvature(sensor::Tuple{Float64,Float64},
				   		  grid::Tuple{Float64,Float64,Float64})
    dist = distance(sensor,(grid[1],grid[2]));
    return grid[3]-curvatureEffect(dist);
end
function correctCurvature(sensor::Tuple{Float64,Float64,Float64},
					   	  grid::Tuple{Float64,Float64,Float64})
    return correctCurvature((sensor[1],sensor[2]),grid);
end

"""
	prismLog(first,second,third)

Auxiliary function used inside sorokin formula (inner product before computing log function)

"""
function prismLog(first::Float64,second::Float64,third::Float64)
    return first + sqrt(first^2 + second^2 + third^2);
end

"""
	prismAtan(first,second,third)

Auxiliary function used inside sorokin formula (inner product + compute atan2)

"""
function prismAtan(first::Float64,second::Float64,third::Float64)
    return atan2(first*sqrt(second^2+third^2+first^2),second*third)
end

end # module
