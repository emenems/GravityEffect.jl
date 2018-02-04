# Unit test for Atmacs gravity effect
function test_atmacsEffect()
	gurl = [pwd()*"\\test\\input\\atmacs_glo.grav"];
	lurl = [pwd()*"\\test\\input\\atmacs_loc.grav"];
	atm1 = atmacsEffect(glofiles=gurl,locfiles=lurl);
	glom = readdlm(gurl[1]);
	locm = readdlm(lurl[1]);
	totalm = -sum(hcat(glom[:,2:end],locm[:,3:end]),2);# - to convert correction to effect
	pressm = locm[:,2]./100;
	timem = DateTime(string.(round.(Int,locm[:,1])),"yyyymmddHH")
	for (i,v) in enumerate(totalm)
		@test atm1[:pressure][i] ≈ pressm[i]
		@test atm1[:effect][i]*1e+9 ≈ v*1e+9
		@test atm1[:datetime][i] == timem[i]
	end
	# Test with interpolation to input time vector
	temp = readdlm(pwd()*"/test/input/atmacs_all_1.grav");
	totalm2 = -sum(temp[:,3:end],2)
	timem2 = DateTime(string.(round.(Int,temp[:,1])),"yyyymmddHH")
	press = DataFrame(datetime=[DateTime(2012,2,29,00,00,00),
								DateTime(2012,2,29,03,00,00)],
						pressure=[temp[1,2]./100,temp[2,2]./100+1])
	atm2 = atmacsEffect(glofiles=[pwd()*"/test/input/atmacs_all_1.grav",
							   pwd()*"/test/input/atmacs_all_2.grav"],
							   pressure=press);
   @test press[:datetime] == atm2[:datetime];
   @test atm2[:effect][1]*1e+9 ≈ totalm2[1]*1e+9
   @test atm2[:effect][2]*1e+9 ≈ totalm2[2]*1e+9-3.0; # -3 for default admittance
end

test_atmacsEffect();
