module Utility

using Printf
"""
Some useful functions
"""

function cutoff_fn(num, min_c, max_c)
    if num < min_c
        return 1.0
    elseif num > max_c
        return 0.0
    else
        t = (num - min_c)/(max_c-min_c)
        return 1.0 - 10.0 * t^3 + 15.0 *  t^4  - 6.0 * t^5
    end
end

function arr2str(a::Array{Int,1})
    st=""
    for i = 1:size(a,1)
        t= @sprintf(" % 2s", a[i])
        st=st*t
    end
    return st

end

function arr2str(a::Array{Float64,1})
    st=""
    for i = 1:size(a,1)
        t= @sprintf(" % 2.10f", a[i])
        st=st*t
    end
    return st

end

#now matches all number arrays!!!
function arr2str(a::AbstractArray{<:Number,2})
    st=""
    for i = 1:size(a,1)
        t= arr2str(a[i,:])
        st=st*t*"\n"
    end
    return st

end


function str_w_spaces(a)
    st=""
    for i = 1:(size(a,1)-1)
        st= st*string(a[i])* "  "
    end
    st= st*string(a[size(a,1)])

    return st

end


function parse_str_ARR_float(sp)

    return map(x->parse(Float64,x),split(sp))

end


function parse_str_ARR_complex(sp)

    
    t = map(x->parse(Float64,x),split(replace(sp, "," => " ")))
    return t[1:2:end] + im*t[2:2:end]

end



end
