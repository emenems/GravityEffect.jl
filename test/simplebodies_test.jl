# bouguerEffect
function bouguerAndCylinder_test()
	thick = 1.234;
	density = 3333.;
	ref = 2.*pi*6.674215e-11*density*thick*1e+9;
	boug_effect = bouguerEffect(thick,density)*1e+9;
	@test ref ≈ boug_effect

	# cylinderEffect
	cylin_effect = cylinderEffect(1.,9e+9,thick,density)*1e+9;
	@test cylin_effect ≈ ref atol=0.01
end

function prismEffect_test()
	thick = 1.234;
	density = 3333.;
	ref = 2.*pi*6.674215e-11*density*thick*1e+9;
	# prismEffect
	prism_effect = prismEffect((0.,0.,1.),(0.,0.,0.),(1.e+7,1.e+7,thick),density)*1e+9;
	@test prism_effect ≈ ref atol=0.1
	prism_effect = prismEffect((0.,0.,0.),(1.,1.,1.),(1.,1.,0.1),2222.)*1e+9;
	@test prism_effect < 0.0
end

# pointEffect
function pointEffect_test()
	point_effect = pointEffect((0.,0.,0.),(0.,1.,0.),1000.)*1e+9;
	@test point_effect == 0.
	point_effect = pointEffect((0.,0.,0.),(1.,1.,10.5),1000.)*1e+9;
	@test point_effect < 0
	ref = 6.674215e-2*3000./4.;# distance 2 m, weight 3000 kg
	point_effect = pointEffect((0.,0.,1.5),(0.,0.,-0.5),3000.)*1e+9;
	@test point_effect ≈ ref atol=0.01
end

function distanceAndCurvature_test()
	# distance2D
	@test distance((0.,0.),(1.,0.)) ≈ 1.
	# distance2D with vector input
	d = distance((1.,2.),[0.,0.],[2.,4.])
	@test d[1] == 1
	@test d[2] ≈ sqrt(1+4)
	# distance3D
	@test distance((0.,0.,0.),(0.,0.,1.)) ≈ 1.
	# prism2point
	@test prism2point((1.,1.,1.),1000.) ≈ 1000.
	# curvatureEffect
	@test curvatureEffect(10000.) ≈ (10/3.57)^2 atol=0.01
	# correctCurvature
	@test correctCurvature((0.,0.),(0.,10000.,10.)) ≈ (10-(10/3.57)^2) atol=0.01
end

bouguerAndCylinder_test();
prismEffect_test();
pointEffect_test();
distanceAndCurvature_test();
