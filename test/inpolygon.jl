@testset "In-polygon" begin
	xv = [0,10,10,0.];
	yv = [0,0,10,10.];
	# Vertex, edge, outside
	x = [0.,10.,3];
	y = [0.,9.9,-0.1];
	o = inpolygon(x,y,xv,yv)
	@test o[1] == true;
	@test o[2] == true;
	@test o[end] == false;
end
