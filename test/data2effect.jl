@testset "Data 2 effect" begin
	resp = DataFrame(layer=[1,2],start=[0.,1.],stop=[1.,2.],
					     total=[4e-9,3e-9],zone1=[3e-9,2e-9],zone2=[1e-9,1e-9])
	sm_in = DataFrame(node1=[0.01,0.2,0.3,0.4],
			node2=[0.00,0.25,0.35,0.45],
			node3=[1.,1.,1.,1.],
			datetime=collect(DateTime(2010,1,1):Dates.Day(1):DateTime(2010,1,4)));
	assign_col = Dict(:node1 => 0.5,:node2 => 1.5, :node3 => 2.5);
	out = sm2effect(sm_in,resp,assign_col)
	@test size(out) == (4,4)
	@test all(isapprox.(out[:node1], sm_in[:node1]*100*resp[:total][1].*1e+9))
	@test all(isapprox.(out[:node2], sm_in[:node2]*100*resp[:total][2].*1e+9))
	@test all(isapprox.(out[:total],out[:node1].+out[:node2]))
	
	
	layer_response = DataFrame(layer = [1,2,3,4], start = [0.0,0.1,0.2,0.3], stop = [0.1,0.2,0.3,0.4],
							total = [4.,50.,61.,70.], zone1 = [3.,40.,50.,60.], zone2 = [1.,10.,11.,10.]);
	layer_agg = aggregate_layer(layer_response,[0.0,0.2,0.4],[0.2,0.4,0.6]);
	# test
	@test layer_agg[:total][1:2] == [4.0+50.0,61.0+70.0]
	@test layer_agg[:zone1][1:2] == [3.0+40.0,50.0+60.0]
	@test layer_agg[:zone2][1:2] == [1.0+10.0,11.0+10.0]
	@test isnan(layer_agg[:total][3])
	@test isnan(layer_agg[:zone1][3])
	@test isnan(layer_agg[:zone2][3])
	@test names(layer_response) == names(layer_agg)
	@test layer_agg[:start] == [0.0,0.2,0.4];
	@test layer_agg[:stop] == [0.2,0.4,0.6];
	@test layer_agg[:layer] == [1,2,3];
end
