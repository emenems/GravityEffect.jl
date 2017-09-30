"""
	atmacsEffect(;glofiles,locfiles,pressure,admittance)
Compute/download Atmacs atmospheric effect on surface gravity

**Input**
* glofiles: see `loadatmacs` function (FileTools Package)
* locfiles:  see `loadatmacs` function (FileTools Package)
* pressure: datafame with :datetime and :pressure columns, where pressure is in hPa
* admittance: admittance used to compute the residual pressure effect in nm/s^2/hPa

**Output**
* DataFrame with DateTime, total effects in m/s^2 and atmospheric pressure in hPa

**Example**
```
# Download data automatically
gurl = ["http://atmacs.bkg.bund.de/data/results/icon/we_icon384_20deg.grav"];
lurl = ["http://atmacs.bkg.bund.de/data/results/iconeu/we_iconeu_70km.grav"];
atm = atmacsEffect(glofiles=gurl,locfiles=lurl)
# Or read data from already downloaded file + resample to local pressure time vector
press = DataFrame(datetime=[DateTime(2012,2,29,00,00,00),
							DateTime(2012,2,29,03,00,00)],
					pressure=[9538.32,9539.32])
e = atmacsEffect(glofiles=["/test/input/atmacs_all_1.grav",
						   "/test/input/atmacs_all_2.grav"],pressure=press)
```
"""
function atmacsEffect(;glofiles::Vector{String}=[""],
					locfiles::Vector{String}=[""],
					pressure::DataFrames.DataFrame=DataFrames.DataFrame(),
					admittance::Float64=-3.0);
	if !isempty(glofiles)
		glo,loc = FileTools.loadatmacs(glofiles=glofiles,locfiles=locfiles);
		if isempty(loc)
			total = glo[:local_newton] .+ glo[:global_newton] .+
						glo[:total_loading];
			press_atmacs = glo[:pressure]./100;
			timevec = glo[:datetime];
		else
			if glo[:datetime] != loc[:datetime]
				glo = ResampleAndFit.interpdf(glo,loc[:datetime])
			end
			total = loc[:local_newton] .+ loc[:regional_newton] .+
					glo[:global_newton] .+ glo[:total_loading];
			press_atmacs = loc[:pressure]./100;
			timevec = loc[:datetime];
		end
		if !isempty(pressure)
			press_atmacs = interpatm(timevec,press_atmacs,pressure[:datetime]);
			total = interpatm(timevec,total,pressure[:datetime]);
			total += (pressure[:pressure] - press_atmacs).*admittance.*1e-9;
			timevec = pressure[:datetime];
		end
		return DataFrames.DataFrame(datetime=timevec,
									effect=total,
									pressure=press_atmacs);
	end
end

"""
Aux function to pass atmacs data for interpolation
"""
function interpatm(x,y,xi);
	return ResampleAndFit.interp1(Dates.value.(x),y,Dates.value.(xi));
end
