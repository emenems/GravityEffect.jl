"""
	eostEffect(filein;timevec)
Compute/download EOST loading effect on surface gravity

**Input**
* filein: file containing the effect or URL to the file
* timevec: time vector used to resample the input effect to required time

**Output**
* DataFrame with DateTime and allo columns of the input file. Warning, EOST loadings are mostly in nm/s^2 (no units conversion)

**Example**
```
# Download data (atmospheric effect) automatically
atm = eostEffect("http://loading.u-strasbg.fr/GGP/atmos/0.10/WE02173h.mog")
# Or read data from already downloaded file (containing Earth rotation effect)
# + resample to given time vector (timein)
timein = @data(collect(DateTime(2012,2,29,00,00,00):Dates.Hour(1):
				DateTime(2012,2,29,03,00,00)));
pol = eostEffect("/test/input/eost_data.rot",timevec=timein);
```
"""
function eostEffect(filein::String;
					timevec::DataFrames.DataArray=DataFrames.@data([]))
	if !isempty(filein)
		eostfile = FileTools.downfile(filein);
		eostdata = FileTools.readggp(eostfile);
		channels = eostchannels(filein);
		if length(names(eostdata)) == length(channels)
			DataFrames.rename!(eostdata,names(eostdata),channels)
		end
		if !isempty(timevec)
			eostdata = ResampleAndFit.interpdf(eostdata,timevec,timecol=channels[1])
		end
		return eostdata;
	end
end

"""
Get channels (fixed order determined by file extension)
"""
function eostchannels(filein::String)
	fileext = split(filein,'.')[end];
	if fileext == "rot"
		return [:datetime,:pol_lod_total,:pol_solid,:pol_ocean,:lod_solid,:lod_load];
	elseif fileext == "oce"
		return [:datetime,:total];
	else #fielex == "hyd" | "inv" | "mog" | "mer"
		return [:datetime,:local,:glob,:total];
	end
end
