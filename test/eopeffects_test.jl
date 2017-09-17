# Unit test for Polar motion and EOP gravity effects
function test_eopeffects()
	eop = eopeffect((0.,90.,0.),filein=pwd()*"/test/input/eop_parameters.c04");
	@test size(eop) == (11,3)
	@test eop[:pol][1]*1e+9 ≈ 0.0
	@test eop[:lod][1]*1e+9 ≈ 0.0
	@test eop[:datetime][2] == DateTime(1962,1,2)
	eop2 = eopeffect((15.,45.,0.),filein=pwd()*"/test/input/eop_parameters.c04",
					amp_factor=0.0);
	@test eop2[:pol][end]*1e+9 ≈ 0.0
	@test eop2[:lod][end]*1e+9 ≈ 0.0
end

test_eopeffects();
