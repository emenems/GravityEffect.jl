"""
	eopEffect(lon,lat,file,amp_factor)
Compute polar motion and length of day effect on surface gravity

**Input**
* sensor: (longitude,latitude,altitude) of the point of computation (deg,deg,m). Altitude will not be used
* filein: either URL to EOP C04 or string with downloaded EOP C04 file
* amp_factor: amplitude factor (default=1.16)

**Output**
* DataFrame with DateTime, Polar motion and Length of Day effects in m/s^2

**Example**
```
# Download data automatically
eop = eopEffect((15.,47.,0.0));
# Or read data from already downloaded file
eop = eopEffect(((15.,47.,0.0),filein="/test/input/eop_parameters.c04",
			amp_factor=1.16);
```
"""
function eopEffect(sensor::Tuple{Float64,Float64,Float64};
			filein::String="http://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now",
			amp_factor::Float64=1.16)
	if !isempty(filein)
		eopfile = FileTools.loadeop(filein);
		# POL: convert to radians
		x = eopfile[:x]./3600.0.*pi./180.0;
		y = eopfile[:y]./3600.0.*pi./180.0;
		# polar motion in m/s^2
		eop = DataFrames.DataFrame(datetime=eopfile[:datetime]);
		eop[:pol] = amp_factor.*poleq.(sensor[1],sensor[2],x,y);
		# LOD: to milisec
		lod = eopfile[:LOD].*1000.0;
		# aux variable
		domega = -0.843994809.*lod./10^12;
		# polar motion in m/s^2
		eop[:lod] = amp_factor.*lodeq.(sensor[1],sensor[2],domega);
		return eop;
	end
end

"""
	poleq(lon,lat,x,y)
Auxiliary funciton to compute polar motion effect
"""
function poleq(lon::Float64,lat::Float64,x::Float64,y::Float64)
	R_const*W_const^2*sind(2.0*lat)*(x*cosd(lon) - y*sind(lon));
end
"""
	lodeq(lon,lat,x,y)
Auxiliary funciton to compute Length of day effect
"""
function lodeq(lon::Float64,lat::Float64,domega::Float64)
	-2*R_const*W_const*cosd(lat)^2*domega;
end
