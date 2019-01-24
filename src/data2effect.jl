"""
	sm2effect(sm,response,assign;mult)
Convert soil moisture to gravity effect using moisture time series and
`layerResponse` output

**Input**
* sm: DataFrame containing soil moisture (in m^3/m^3) time series
* response: output of `layer2response`
* assign: assign sm names to depths
* mult: multiply the effect by given factor, e.g. to get results in nm/s^2. Response expected to be in m/s^2, thus mult=1e+9 to get nm/s^2

**Output**
* DataFrame containing results in nm/s^2 (or units specified by `mult` parameter)

**Example**
```
resp = DataFrame(layer=[1,2],start=[0.,1.],stop=[1.,2.],
				     total=[4e-9,3e-9],zone1=[3e-9,2e-9],zone2=[1e-9,1e-9])
sm_in = DataFrame(node1=[0.1,0.2,0.3,0.4],
		node2=[0.15,0.25,0.35,0.45],
		datetime=collect(DateTime(2010,1,1):Dates.Day(1):DateTime(2010,1,4)));
assign_col = Dict(:node1 => 0.5,:node2 => 1.5);
sm2effect(sm_in,resp,assign_col)
```
"""
function sm2effect(sm::DataFrame,response::DataFrame,
				   assign::Dict{Symbol,Float64};mult::Float64=1e+9)::DataFrame
	out = DataFrame(datetime=sm[:datetime],
					total=zeros(Float64,length(sm[:datetime])));
	for i in keys(assign)
		r = finddepth(assign[i],response[:start],response[:stop])
		if !isempty(r)
			# multiply by 100 to get Vol% & convert to nm/s^2 (use only one SM
			# time series even if more found => r[1])
			out[i] = response[:total][r[1]].*sm[i].*100.0.*mult;
			out[:total] += out[i];
		end
	end
	return out
end

"""
Auxiliary function to find index (in response) corresponding to given depth
"""
function finddepth(currdepth::Float64,
			startdepths::Vector{Float64},stopdepths::Vector{Float64})
	r = map(x->x<=currdepth,startdepths) .& map(x->x>currdepth,stopdepths)
	return findall(x->x.==true,r) # return index
end

"""
Function to aggregate high-resolution unit gravity response to input depth range

**Input**
* gin: result of the `layerResponse` function
* dstart: starting depth (same units as in `gin`)
* dstop: end depth (same units as in `gin`)

**Output**
* DataFrame in the same format as `layerResponse` output (=gin input)

**Example**
```
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
```
"""
function aggregate_layer(gin::DataFrame,dstart::Vector{Float64},dstop::Vector{Float64})::DataFrame
	# declare output
	df = DataFrame(zeros(length(dstart),ncol(gin)).+NaN, names(gin))
	df[:layer] = collect(1:length(dstart))
	df[:start] = dstart;
	df[:stop] = dstop;
    for i in 1:length(dstart)
        r = map(x->x>=dstart[i],gin[:start]) .& map(x->x<=dstop[i],gin[:stop]);
        if any(r)
            for z in 4:ncol(gin)
                df[z][i] = sum(gin[z][r])
            end
        end
    end
    return df
end
