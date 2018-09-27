# layerResponse
function layerResponse_test1()
	sensor=Dict(:x => 100., :y => 150., :z => 100., # coordinates (z=altitude)
				:sensHeight => 1.0) # will be added to altitude (:z)
	layers=Dict(:start => [0.0, 1.0], # starting depth of the soil layer
				:stop  => [1.0, 2.0]) # end of the soil layer
	dem_in = pwd()*"/test/input/dem_data.asc";
	zones =Dict(:dem   => [dem_in,dem_in], # input DEMs
				:radius=> [50.,200.], # integration radius
				:resolution => [0.1,0.2], # resolution of the zone (set to negative value to keep original DEM resolution)
				:interpAltitude => [true,false]); # interpolate the altitude of the sensor (== sensor[:z])
	exclude = Dict(:cylinder=>Dict(:radius=>1.,:start=>0.5,:stop=>2.),
				   :prism=>Dict(:x=>[100.,106.],:y=>[155.,156.],:dx=>[2.,3.],:dy=>[1.5,1.8],
				   				:start=>[0.,1.],:stop=>[1.,2.]),
				   :polygon=>Dict(:file=>pwd()*"/test/input/exclusion_polygon.txt",
				   				  :start=>0.,:stop=>2.));
	outdata0 = layerResponse(sensor,layers,zones)
	outdata1 = layerResponse(sensor,layers,zones,exclude=exclude,
				outfile=pwd()*"/test/output/layerResponse_results.txt");

	# compute "should" value
	mainpart = [prismEffect((sensor[:x],sensor[:y],sensor[:z]+sensor[:sensHeight]),
						(sensor[:x],sensor[:y],sensor[:z]-layers[:start][1]),
						(150.0+zones[:resolution][2],250.0+zones[:resolution][2],1.),10.),
			    prismEffect((sensor[:x],sensor[:y],sensor[:z]+sensor[:sensHeight]),
			 			(sensor[:x],sensor[:y],sensor[:z]-layers[:start][2]),
						(150.0+zones[:resolution][2],250.0+zones[:resolution][2],1.),10.)]
	zone1 = [cylinderEffect(sensor[:sensHeight]+layers[:start][1],zones[:radius][1],
						layers[:stop][1]-layers[:start][1],10.),
			cylinderEffect(sensor[:sensHeight]+layers[:start][2],zones[:radius][1],
		   		 		layers[:stop][2]-layers[:start][2],10.)];
	exclude1 = [cylinderEffect(sensor[:sensHeight]+exclude[:cylinder][:start],
						exclude[:cylinder][:radius],
						layers[:start][2]-exclude[:cylinder][:start],10.),
			    cylinderEffect(sensor[:sensHeight]+layers[:start][2],
						exclude[:cylinder][:radius],
						layers[:stop][2]-layers[:start][2],10.)];
	exclude2 = [prismEffect((sensor[:x],sensor[:y],sensor[:z]+sensor[:sensHeight]),
						(exclude[:prism][:x][1],exclude[:prism][:y][1],sensor[:z]-layers[:start][1]),
						(exclude[:prism][:dx][1],exclude[:prism][:dy][1],layers[:stop][1]-layers[:start][1]),10.0),
				prismEffect((sensor[:x],sensor[:y],sensor[:z]+sensor[:sensHeight]),
						(exclude[:prism][:x][2],exclude[:prism][:y][2],sensor[:z]-layers[:start][2]),
						(exclude[:prism][:dx][2],exclude[:prism][:dy][2],layers[:stop][2]-layers[:start][2]),10.0)];
	# polysong (see the input file)
	exclude3 = [prismEffect((sensor[:x],sensor[:y],sensor[:z]+sensor[:sensHeight]),
						(25.0+20.0,225.0,sensor[:z]-layers[:start][1]),
						(40.0,70.0,layers[:stop][1]-layers[:start][1]),10.)
				prismEffect((sensor[:x],sensor[:y],sensor[:z]+sensor[:sensHeight]),
						(25.0+20.,225.0,sensor[:z]-layers[:start][2]),
						(40.0,70.0,layers[:stop][2]-layers[:start][2]),10.0)]
	compdata = (mainpart-exclude1-exclude2-exclude3).*1e+9

	# Compare results without exclusion
	@test isapprox(outdata0[:total].*1e+9,mainpart.*1e+9,atol=1e-4)
	@test isapprox(outdata0[:zone1].*1e+9,zone1.*1e+9,atol=1e-4)

	# Compare results with exclusion bodies
	@test isapprox(outdata1[:total].*1e+9,compdata,atol=1e-4)
end

function layerResponse_test_aux1()
	dem_in = pwd()*"/test/input/dem_data.asc";
	dem = GravityEffect.response_load_dem(dem_in)
	@test dem[:y] == collect(25.:50:25.0+5*50)
	@test dem[:x] == collect(25.:50:25.0+3*50)
	@test dem[:height] == zeros(Float64,(length(dem[:y]),length(dem[:x]))).+100.0
	@test GravityEffect.interpheight(26.,76.,dem) == 100.
	@test isnan(GravityEffect.interpheight(24.99,76.,dem))

 	depth_start,depth_stop = 2.,3.
	comp,ubound,lbound = GravityEffect.iswithin(2.5,3.1,depth_start,depth_stop)
	@test comp
	@test lbound==3.0
	@test ubound==2.5
	comp,ubound,lbound = GravityEffect.iswithin(0.5,2.,depth_start,depth_stop)
	@test !comp
	sensor = Dict(:x=>100.,:y=>200.,:zUse=>101.5)
	prism = Dict(:x=>100.,:y=>210.,:start=>2.,:stop=>3.,:dx=>1.,:dy=>2.)
	prismcopy = deepcopy(prism);prismcopy[:density] = 1111.;prismcopy[:z] = NaN;
	dem[:x],dem[:y] = ResampleAndFit.meshgrid(dem[:x],dem[:y])
	GravityEffect.interpdata!(prismcopy,dem[:x],dem[:y],dem[:height],:z);
	@test prismcopy[:z] == 100.

	# Excluding prism (full layer)
	tp = GravityEffect.excludeprism(sensor,prismcopy,depth_start,depth_stop).*1e+9
	dg = prismEffect((0.,0.,1.5),(0.,10.,-depth_start),(1.,2.,1.),1111.)*1e+9;
	@test tp ≈ dg
	# For prism outside of layer
	tp = GravityEffect.excludeprism(sensor,prismcopy,4.,5.).*1e+9
	@test tp ≈ 0.0

	# Excluding cylinder (half of the layer )
	cylinder = Dict(:start=>2.5,:stop=>3.5, :radius=>1.,
					:x=>sensor[:x],:y=>sensor[:y])
	GravityEffect.interpdata!(cylinder,dem[:x],dem[:y],dem[:height],:z);cylinder[:density] = 1000.
	tc = GravityEffect.excludecylinder(sensor,cylinder,depth_start,depth_stop).*1e+9
	dg = cylinderEffect(1.5+cylinder[:start],cylinder[:radius],0.5,1000.)*1e+9;
	@test tc ≈ dg
end
function layerResponse_test_aux2()
	dem_in = pwd()*"/test/input/dem_data.asc";
	dem = GravityEffect.response_load_dem(dem_in)
	dem[:x],dem[:y] = ResampleAndFit.meshgrid(dem[:x],dem[:y]);
	polyg = Dict(:file=>pwd()*"/test/input/exclusion_polygon.txt",
				 :start=>0.,:stop=>1.);
	density0 = dem[:x].*0.0 .+ 10.0;
	density = copy(density0);
	GravityEffect.excludepolygon!(polyg,dem,0.,1.,density)
	@test sum(density) == 10.0*(size(dem[:x],1)*size(dem[:y],2))-10.0
	@test density[5,1] == 0.0
	density = copy(density0);
	GravityEffect.excludepolygon!(polyg,dem,1.,2.,density)
	@test sum(density) == 10.0*(size(dem[:x],1)*size(dem[:y],2))
end

layerResponse_test_aux1();
layerResponse_test_aux2();
layerResponse_test1();
