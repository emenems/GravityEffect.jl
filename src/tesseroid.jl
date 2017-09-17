"""
Calculate the gravitational effect of a tesseroid of given position.
Calculation is based on Heck and Seitz paper "A comparison of tesseroid,
prism and point-mass approaches from mass reduction in gravity field
modelling", 2007, J Geodesy, 81:121-136, (doi:10.1007/s00190-006-0094-0).
It is a spherical approximation (Taylor expansion) not suitable for near
zone (see paper)! Use higher resolution for local zone.
The code has not been optimized yet!

**Input**
* sensor: (longitude,latitude,altitude) of the point of computation (deg,deg,m)
* tesser: (longitude,latitude,altitude) of the tesseroid (deg,deg,m). altitude should point to LOWER boundary, while long and lat to center of the tesseroid.
* resol: (dLon,dLat,heihgt) resolution of the tesseroid (deg,deg,m)
* density: differential density of the tesseroid in kg/m^3 (by default 1.0)
* radius: radius of the reference sphere (in m, by default Module constatn)

**Output**
* resulting gravity effect in m/s^2

**Example**
```
dg = tesseroid((0.,0.,10.),(1.,1.,0.),(0.01,0.01,1.),1000.,6371000.)
```
"""
function tesseroid(sensor::Tuple{Float64,Float64,Float64},
				   tesser::Tuple{Float64,Float64,Float64},
				   resol::Tuple{Float64,Float64,Float64},
				   density::Float64=1.0,
				   radius::Float64=R_const)
	# Prepare computation point=sensor
	lon = deg2rad(sensor[1]);
	lat = deg2rad(sensor[2]);
	sinlat = sin(lat);
	coslat = cos(lat);
	r = radius + sensor[3];
	r2 = r^2;

	# Running integration point
	lon0 = deg2rad(tesser[1]);
	lat0 = deg2rad(tesser[2]);
	coslat0 = cos(lat0);
	sinlat0 = sin(lat0);
	sinlat02 = sinlat0^2;
	r0 = radius + tesser[3] + resol[3]/2;
	r02 = r0^2;
	r03 = r0^3;
	delta_r = resol[3];
	delta_lonT = deg2rad(resol[1]);
	delta_latT = deg2rad(resol[2]);
	delta_lon = lon0-lon;
	cosdelta_lon = cos(delta_lon);
	sindelta_lon = sin(delta_lon);
	coslat0xcosdelta_lon = coslat0*cosdelta_lon;
	coslatxsinlat0xcosdelta_lon = coslat*sinlat0*cosdelta_lon;
	sinlatxsinlat0 = sinlat*sinlat0;
	sinlatxcoslat0 = sinlat*coslat0;

	# Calc auxiliary variables
	psi0 = acos(sinlatxsinlat0 + coslat*coslat0*cos(lon0 - lon));
	cospsi0 = cos(psi0);
	l0 = sqrt(r2 + r02 - 2*r*r0*cospsi0);
	l02 = l0^2;
	l03 = l0^3;
	l04 = l02*l02;
	r0dl03 = (r0/l0)^3;

	# Calc integral kernels
	L000 = (r02*(r - r0*cospsi0)*coslat0)/l03;

	L200 = ((r*coslat0)/l03)*(2-((3*r0)/l02)*(5*r0 - (2*r + 3*r0*cospsi0)*cospsi0) +
	        ((15*r03)/l04)*sin(psi0)^2*(r0 - r*cospsi0));

	L020 = r0dl03*coslat*(1 - 2*sinlat02)*cosdelta_lon +
	        (r02/l0^5)*(-r*(r2 + r02)*coslat0 +
	        r0*sinlat*(-r*r0*(sinlatxcoslat0 - coslatxsinlat0xcosdelta_lon) +
	        sinlat0*coslat0*(2*r2 + 4*r02 - 3*r*r0*sinlatxsinlat0)) +
	        r02*coslat*cosdelta_lon*(1 - 2*sinlat02)*
	        (r0 + r*coslat*coslat0xcosdelta_lon) +
	        r*r02*coslat*sinlat0*coslat0xcosdelta_lon*
	        (3*sinlatxcoslat0 - 4*coslatxsinlat0xcosdelta_lon)) +
	        ((5*r*r03)/l0^7)*(-r*(r2 + r02)*sinlat0 +
	        r02*coslat*sinlat0*coslat0xcosdelta_lon*
	        (r0 + r*coslat*coslat0xcosdelta_lon) +
	        r0*sinlat*(2*r2 - r02 - r*r0*cospsi0 + sinlat02*
	        (r2 + 2*r02 - r*r0*sinlatxsinlat0)))*
	        (sinlatxcoslat0 - coslatxsinlat0xcosdelta_lon);

	L002 = r0dl03*coslat*coslat0^2*
	        (cosdelta_lon - ((3*r)/l02)*(2*r0*coslat*coslat0*sindelta_lon^2 +
	        (r - r0*cospsi0)*cosdelta_lon) +
	        ((15*r2*r0)/l04)*coslat*coslat0*(r - r0*cospsi0)*sindelta_lon^2);

	# Final calc
	return G_const*density*delta_r*delta_latT*delta_lonT*(L000 + (1/24)*(L200*delta_r^2 +
	 		L020*delta_latT^2 + L002*delta_lonT^2));
end

function deg2rad(deg::Float64)
	return deg*pi/180;
end
