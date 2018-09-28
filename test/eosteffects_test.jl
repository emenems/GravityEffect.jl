# Unit test for EOST Loading effects
function test_eostEffects()
	# Load without interpolation
	eost1 = eostEffect(joinpath(dirname(@__DIR__),"test","input","eost_data.rot"));
	eost1m = DelimitedFiles.readdlm(joinpath(dirname(@__DIR__),"test","input","eost_data.rot"),skipstart=14)[1:end-1,:];
	@test names(eost1) == [:datetime,:pol_lod_total,:pol_solid,
							:pol_ocean,:lod_solid,:lod_load];
	@test eost1[:datetime] == DateTime.(string.(eost1m[:,1]),"yyyymmdd",locale="english")
	for i in 3:size(eost1m,2)
		@test sum(eost1[i-1]) ≈ sum(eost1m[:,i])
	end

	# Load + interpolation
	timein = [DateTime(1980,01,02)];
	eost2 = eostEffect(joinpath(dirname(@__DIR__),"test","input","eost_data.rot"),timevec=timein);
	@test size(eost2) == (1,6)
	@test eost2[:datetime] == timein
	for i in 3:size(eost1m,2)
		@test eost2[i-1][1] ≈ eost1m[2,i]
	end
end

test_eostEffects();
