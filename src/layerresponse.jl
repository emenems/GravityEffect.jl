"""
layerResponse(sensor,layers,zones;exclude,nanheight,outfile,def_density)
Compute gravity response for soil layers copying the terrain

**Input**
* sensor: Dictionary containing information about the position of the sensor (see example below)
* layers: Dictionary with soil layer specifications (=>:start & :stop depths), i.e. soil layers copying the terrain
* zones: Dictionary setting the computation zones separated by integration radii (see example below)
* nanheight: return values even when DEM either contains NaNs or the integration radius exceeds the DEM range (`true` by default)
* exclude: exclude bodies from computation (first zone only!): prism, cylinder or all points inside a polygon.
* outfile: output file (no output file by default)
* def_density: default density in kg/m^3 (10 => 1Vol% == 0.01 m^3/m^3)

**Important notes**
When using multiple bodies in the `exclude` parameter, pay attention to possible overlapping.
Overlapping will lead to incorrect subtraction of the body.
Set exclusion bodies (prism and cylinder) only for the first zone!
The exclusion depths do not need to coincide with those in `layers` (e.g. part of the body will be subtracted from soil layer)
To convert Shapefile to supported text file, use the Octave code in /src/auxiliary/shp2txt.m
Use only "simple" polygon, i.e. not polygon inside a polygon.
Density of 10 kg/m^3 corresponds to Volumetric content of 1%.
The effect of the Earth's curvature on height (input DEM) will be automatically corrected.
Results were (successfully) compared to https://www.hydrol-earth-syst-sci.net/21/3167/2017/

**Output**
* DataFrame containing results for each layer and each zone in m/s^2

**Example**
```
sensor=Dict(:x => 100., :y => 150., :z => 100., # coordinates (z=altitude)
		:sensHeight => 1.0) # will be added to altitude (:z)
layers=Dict(:start => [0.0, 1.0], # starting depth of the soil layer (relative to surface)
		:stop  => [1.0, 2.0]) # end of the soil layer (relative to surface)
dem_in = pwd()*"/test/input/dem_data.asc";
zones =Dict(:dem   => [dem_in,dem_in], # input DEMs (identical for both zones in this case )
		:radius=> [50.,200.], # integration radius
		:resolution => [0.1,0.2], # resolution of the zone (linear interpolation, set to negative value to keep original DEM resolution)
		:interpAltitude => [true,false]); # linearly interpolate the altitude of the sensor (== sensor[:z])
exclude = Dict(:cylinder=>Dict(:radius=>1.,:start=>0.5,:stop=>2.), # one cylinder. Must be below sensor (no :x,:y), start and stop relative to surface
		   :prism=>Dict(:x=>[100.,106.],:y=>[155.,156.], # (absolute) coordinates of the prism
						:dx=>[2.,3.],:dy=>[1.5,1.8], # prism size
						:start=>[0.,1.],:stop=>[1.,2.]), # relative starting and end depths
		   :polygon=>Dict(:file=>pwd()*"/test/input/exclusion_polygon.txt", # exclusion polygon file
						  :start=>0.,:stop=>2.)) # relative starting and end depths
outputfile = pwd()*"/test/output/layerResponse.txt";
outdata = layerResponse(sensor,layers,zones,
					exclude=exclude,nanheight=true,
					outfile=outputfile,def_density=10.);
```
"""
function layerResponse(sensor::Dict{Symbol,Float64},
					 layers::Dict{Symbol,Array{Float64,1}},
					 zones::Dict{Symbol,Any};
					 exclude::Dict{Symbol,}=Dict(:nothing=>""),
					 nanheight::Bool=true,
					 outfile::String="",
					 def_density::Float64=10.)::DataFrame
	 # declare output variable
	 outdata = DataFrames.DataFrame(layer = collect(1:1:length(layers[:start])),
	 					 start = layers[:start],stop = layers[:stop],
						 total = zeros(Float64,length(layers[:start])));
	 ## run loop for all zones
	 for i in 1:length(zones[:dem])
		 # append zoneX to output dataframe (fill with dummies)
		 zoneX = Symbol("zone"*string(i)); # name of the current zone
		 outdata[zoneX] = zeros(Float64,length(layers[:start]));
		 # get/prepare current DEM
		 dem = response_prepare_dem(zones[:dem][i],zones[:radius][i],
									zones[:resolution][i],sensor);
		 # interpolate sensor altitude if needed & add sensor height
		 sensor[:zUse] = zones[:interpAltitude][i] ? interpheight(sensor[:x],sensor[:y],dem) : sensor[:z];
		 sensor[:zUse] += sensor[:sensHeight];
		 # get height of exclusion body (first zone only)
		 if haskey(exclude,:cylinder) && i == 1
			 cylinders = deepcopy(exclude[:cylinder]);
			 # use for interpolation
			 cylinders[:x] = length(cylinders[:stop])!=1 ? zeros(length(cylinders[:stop])).+sensor[:x] : sensor[:x];
			 cylinders[:y] = length(cylinders[:stop])!=1 ? zeros(length(cylinders[:stop])).+sensor[:y] : sensor[:y];
			 interpdata!(cylinders,dem[:x],dem[:y],dem[:height],:z)
		 end
		 if haskey(exclude,:prism) && i == 1
			 prisms = deepcopy(exclude[:prism]);
			 interpdata!(prisms,dem[:x],dem[:y],dem[:height],:z)
		 end
		 # density matrix used to include only grid cells within current zone
		 density = zeros(Float64,size(dem[:height])) .+ def_density;
		 density[dem[:distance].>zones[:radius][i]] = 0.; # all points above intergration radius
		 if i > 1 # remove adjacent, i.e. the previous zone
			 density[dem[:distance].<=zones[:radius][i-1]] = 0.;
		 end
		 if nanheight
			 density[isnan.(dem[:height])] .= 0.0;
			 dem[:height][isnan.(dem[:height])] .= 9e+30;
		 end
		 ## loop for each soil layer adding all grid cells
		 for j in 1:length(layers[:start])
			 density_use = deepcopy(density);
			 if haskey(exclude,:polygon)
				 excludepolygon!(exclude[:polygon],dem,
				 				layers[:start][j],layers[:stop][j],
								density_use);
			 end
			 # first compute effect including all points
			 outdata[zoneX][j] = response_grid_cells(sensor,dem,density_use,
			 								layers[:start][j],layers[:stop][j]);
			 # exclude (subtract) requested bodies (will be applied to the first
   			 # computation zone (i==1) only)
			 if i == 1
				 if haskey(exclude,:cylinder)
					 interpdata!(cylinders,dem[:x],dem[:y],density_use,:density);
					 outdata[zoneX][j] -= excludecylinder(sensor,cylinders,
				 								layers[:start][j],layers[:stop][j]);
				 end
				 if haskey(exclude,:prism)
					 interpdata!(prisms,dem[:x],dem[:y],density_use,:density);
					 outdata[zoneX][j] -= excludeprism(sensor,prisms,
 				 								layers[:start][j],layers[:stop][j]);
				 end
			 end
		 end
		 # add to total effect
		 outdata[:total] += outdata[zoneX];
	 end
	 if !isempty(outfile)
		 FileTools.write_layerResponse(sensor,layers,zones,exclude,nanheight,
		 						outfile,def_density,outdata)
	 end
	 return outdata
end

"""
auxiliary function to prepare DEM for computation: load + resample to computation grid
"""
function response_prepare_dem(file_in::String,radius_in::Float64,resol_in::Float64,
							  sensor::Dict{Symbol,Float64})::Dict{Symbol,Any}
	dem = response_load_dem(file_in);
	# resample DEM
	if resol_in > 0
		x,y = ResampleAndFit.meshgrid(create_dem_vector(sensor[:x],resol_in,radius_in),
									  create_dem_vector(sensor[:y],resol_in,radius_in));
		demi = Dict(:x=>x,:y=>y,:height=>ResampleAndFit.interp2(dem[:x],dem[:y],dem[:height],
												  x,y));
	else
		x,y = ResampleAndFit.meshgrid(dem[:x],dem[:y]);
		demi = Dict(:x=>x,:y=>y,:height=>dem[:height]);
	end
	# compute distance to the sensor
	demi[:distance] = distance((sensor[:x],sensor[:y]),demi[:x],demi[:y]);
	# correct for Earth's curvature
	demi[:height] -= curvatureEffect.(demi[:distance]);
	return demi
end

"""
auxiliary function to load DEMs
"""
function response_load_dem(file_in::String)::Dict{Symbol,Any}
	if split(file_in,".")[end] == "asc"
		return FileTools.loadascii(file_in);
	end
end

"""
auxiliary function to create new DEM vector (for interpolation)
"""
function create_dem_vector(x::Float64,dx::Float64,r::Float64)::Vector{Float64}
	# make sure the extantion is bigger than the radius (=center+radius+2*resolution)
	return collect((x-r-dx*2):dx:(x+r+dx*2));
end

"""
auxiliary function to interpolate height of a point
"""
function interpheight(x::Float64,y::Float64,dem::Dict{Symbol,Any})::Float64
	ResampleAndFit.interp2(dem[:x],dem[:y],dem[:height],x,y)
end

"""
auxiliary function to produce a sum of all gravity effects for current (=one)
layer
"""
function response_grid_cells(sensor::Dict{Symbol,Float64},
						    dem::Dict{Symbol,Any},density::Matrix{Float64},
							start_depth::Float64,end_depth::Float64)::Float64
	# get constatnt resolution of the DEM & coordinates of the point of computation
	resol = (dem[:x][1,2]-dem[:x][1,1], # dx
			 dem[:y][2,1]-dem[:y][1,1], # dy
		     end_depth-start_depth); # dz
   	xyz = (sensor[:x],sensor[:y],sensor[:zUse]);
	# declare output
	out::Float64 = 0.0;
	# run loop for all grid cells
	for (i,v) in enumerate(density[:])
		if !isapprox(v,0.0)
			out += prismEffect(xyz,
						  (dem[:x][i],dem[:y][i],dem[:height][i]-start_depth), # i-th prism coordinates
						  resol,v); # size of the prism + its density
		end
	end
	return out
end

"""
Auxiliary function to compute the gravity effect of a cylinder inside a layer
Input Dictionary can contain multiple cylinders (:radius = vector)
"""
function excludecylinder(sensor::Dict{Symbol,Float64},cylinder::Dict{Symbol,},
						start_depth::Float64,end_depth::Float64)::Float64
	out::Float64=0.0;
	for i in 1:length(cylinder[:radius])
		comp,ubound,lbound = iswithin(cylinder[:start][i],cylinder[:stop][i],start_depth,end_depth)
		if comp
			out+=cylinderEffect((sensor[:zUse]-cylinder[:z][i])+ubound, # starting depth with respect to sensor
							cylinder[:radius][i],lbound-ubound,cylinder[:density][i]) #radius, thickness and density
		end
	end
	return out
end

"""
Auxiliary function to compute the gravity effect of a prism inside a layer
Input Dictionary can contain multiple prisms (:xyz & :dxdydz = vectors)
"""
function excludeprism(sensor::Dict{Symbol,Float64},prism::Dict{Symbol,},
					    start_depth::Float64,end_depth::Float64)::Float64
	out::Float64=0.0;
	for i in 1:length(prism[:x])
		comp,ubound,lbound = iswithin(prism[:start][i],prism[:stop][i],start_depth,end_depth)
		if comp
			out += prismEffect((sensor[:x],sensor[:y],sensor[:zUse]), # computation point
						(prism[:x][i],prism[:y][i],prism[:z][i]-ubound), # prism position
						(prism[:dx][i],prism[:dy][i],lbound-ubound),prism[:density][i]); # prism size
		end
	end
	return out
end
"""
Auxiliary function to get density/height in the layer of the exclusion prism/cylinder
"""
function interpdata!(bodyin::Dict{Symbol,},
				x::Matrix{Float64},y::Matrix{Float64},d::Matrix{Float64},
				col::Symbol)
	temp = [];
	for i in 1:length(bodyin[:x])
		push!(temp,ResampleAndFit.interp2(x,y,d,bodyin[:x][i],bodyin[:y][i]))
	end
	bodyin[col] = length(bodyin[:x])==1 ? temp[1] : temp;
end

function excludepolygon!(polyg::Dict{Symbol,},dem::Dict{Symbol,Any},
						start_depth::Float64,end_depth::Float64,
						density_use::Matrix{Float64})
	for (i,v) in enumerate((length(polyg[:start])==1 ? [polyg[:file]] : polyg[:file]))
		if iswithin(polyg[:start][i],polyg[:stop][i],start_depth,end_depth)[1]
			polygxy = readpolygon(v);
			# one file can contain multiple polygons separated by NaNs
			r = any(isnan.(polygxy[:x])) ? find(isnan,polygxy[:x]) : [length(polygxy[:x])+1];
			s = 1;
			for j in r
				density_use[inpolygon(dem[:x],dem[:y],polygxy[:x][s:j-1],polygxy[:y][s:j-1])] = 0.;
				s = j+1;
			end
		end
	end
end
function readpolygon(filein)
	temp = readdlm(filein,comment_char='%');
	DataFrame(x = temp[:,1],y = temp[:,2]);
end
"""
check if exclusion body (cylinder or prism) reach up to the current soil layer
and return upper and lower boundaty of the exclusion body
"""
function iswithin(bodystart,bodystop,start_depth,end_depth)
	upperbound = bodystart<=start_depth ? start_depth : bodystart;
	lowerbound = bodystop>=end_depth ? end_depth : bodystop;
	return lowerbound>upperbound,upperbound,lowerbound
end
