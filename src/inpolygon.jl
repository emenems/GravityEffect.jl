"""
    inpolygon(x,y,xv,yv)

Check if point is inside a polygon exploiting Octave [function](
https://octave.sourceforge.io/octave/function/inpolygon.html)

**Input**
* x: scalar, vector or matrix of x coordinate to be identified
* y: scalar, vector or matrix of y coordinate to be identified
* xy: vector of the polygon x coordinates
* xy: vector of the polygon y coordinates

**Example**

```
id = inpolygon(0.2,0.3,[0.,1.,1.,0.,0.],[0.,0.,1.,1.,0.])

idvec = inpolygon([0.2,0.0],[0.3,0.0],[0.,1.,1.,0.,0.],[0.,0.,1.,1.,0.])

```

"""

function inpolygon{T<:Real}(x::T, y::T, xv::Vector{T}, yv::Vector{T})
    npol = j = length(xv);
    #inside = on = Vector{Bool}(length(x));
	inside = on = false;
    for i = 1:npol
        delta_xv = xv[j] - xv[i];
        delta_yv = yv[j] - yv[i];
        # distance = [distance from (x,y) to edge] * length(edge)
        distance = delta_xv * (y - yv[i]) - (x - xv[i]) * delta_yv;
		# is y between the y-values of edge i,j AND (x,y) on the left of the edge?
    	idx1 = ((((yv[i] <= y) & (y < yv[j])) | ((yv[j] <= y) & (y < yv[i])))
				& (0. < distance*delta_yv));
		if idx1
    		inside = !inside;
		end
    	# Check if (x,y) are actually on the boundary of the polygon.
    	idx2 = ((((yv[i] <= y) & (y <= yv[j])) | ((yv[j] <= y) & (y <= yv[i])))
              & (((xv[i] <= x) & (x <= xv[j])) | ((xv[j] <= x) & (x <= xv[i])))
    	 	  & ((0. == distance) | (delta_xv == 0.)));
		if idx2
    		on = true;
		end
	    j = i;
  	end
	# Matlab/Octave definition includes both in polygon and on polygon points.
  	return (inside .| on);
end
# Version for vector input
function inpolygon{T<:Real}(x::Vector{T}, y::Vector{T}, xv::Vector{T}, yv::Vector{T})
    id = Vector{Bool}(length(x));
    for i = 1:length(x)
        id[i] = inpolygon(x[i],y[i],xv,yv);
    end
    return id
end
# Version for matrix input
function inpolygon{T<:Real}(x::Matrix{T}, y::Matrix{T}, xv::Vector{T}, yv::Vector{T})
    id = inpolygon(x[:],y[:],xv,yv);
    return reshape(id,size(x))
end
