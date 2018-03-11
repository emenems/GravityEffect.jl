"""
	sm2effect(sm,response,assign)
Convert soil moisture to gravity effect using moisture time series and
`layerResponse` output

**Input**
* sm: DataFrame containing soil moisture (in m^3/m^3) time series
* response: output of `layer2response`
* assign: assign sm names to depths 

**Output**
* DataFrame containing results in nm/s^2

**Example**
```
resp = DataFrame(layer=[1,2],start=[0.,1.],stop=[1.,2.],
				     total=[4e-9,3e-9],zone1=[3e-9,2e-9],zone2=[1e-9,1e-9])
sm_in = DataFrame(node1=[0.1,0.2,0.3,0.4],
		node2=[0.15,0.25,0.35,0.45],
		datetime=collect(DateTime(2010,1,1):Dates.Day(1):DateTime(2010,1,4)));
assign_col = Dict(:node1 => 0.5,:node2 => 1.5);
```
"""
function sm2effect(sm::DataFrame,response::DataFrame,
				   assign::Dict{Symbol,Float64})::DataFrame
	out = DataFrame(datetime=sm[:datetime],
					total=zeros(Float64,length(sm[:datetime])));
	for i in keys(assign)
		r = finddepth(assign[i],response[:start],response[:stop])
		if !isempty(r)
			# multiply response by SM and convert to nm/s^2 (use only one SM
			# time series even if more found => r[1])
			out[i] = response[:total][r[1]].*sm[i].*1e+9;
			out[:total] += out[i];
		end
	end
	return out
end

"""
Auxiliary function to find index (in response) corresponding to given depth
"""
function finddepth(currdepth::Float64,
			startdepths::DataArray{Float64},stopdepths::DataArray{Float64})
	r = map(x->x<=currdepth,startdepths) .& map(x->x>currdepth,stopdepths)
	return find(x->x.==true,r) # return index
end