function get_ll_c_func(ξ, _ll::Function)::Function
    return (x) -> _ll(x, ξ)
end

function get_ll_ξ_func(c::Vector{Vector{Float64}}, _ll::Vector{Function})
    _f = (x) -> [_ll[i](c[i], x) for i in eachindex(c)]
    _fsum = x -> sum(_f(x))
    return _f, _fsum
end

function get_ll_func_c_ξ(dataobs::Matrix{Float64}, dt)
    nind = size(dataobs)[2]
    ll = Vector{Function}(undef, nind)
    for i in 1:nind
        ll[i] = (log_c, log_ξ) -> kalman_filter(dataobs[:, i], log_ξ, log_c, dt)
    end
    return ll
end
