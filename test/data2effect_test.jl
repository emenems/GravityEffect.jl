function sm2effect_test()
	resp = DataFrame(layer=[1,2],start=[0.,1.],stop=[1.,2.],
					     total=[4e-9,3e-9],zone1=[3e-9,2e-9],zone2=[1e-9,1e-9])
	sm_in = DataFrame(node1=[0.1,0.2,0.3,0.4],
			node2=[0.15,0.25,0.35,0.45],
			node3=[1.,1.,1.,1.],
			datetime=collect(DateTime(2010,1,1):Dates.Day(1):DateTime(2010,1,4)));
	assign_col = Dict(:node1 => 0.5,:node2 => 1.5, :node3 => 2.5);
	out = sm2effect(sm_in,resp,assign_col)
	@test size(out) == (4,4)
	@test all(isapprox.(out[:node1], sm_in[:node1]*resp[:total][1].*1e+9))
	@test all(isapprox.(out[:node2], sm_in[:node2]*resp[:total][2].*1e+9))
	@test all(isapprox.(out[:total],out[:node1].+out[:node2]))
end

sm2effect_test();
