function globalset(;start=1e-3u"kpc", radius=16.0u"kpc", datapoints::Int=100, mode::Symbol=:linear)
    params = SetupParams()

    params.mode=mode
    params.datapoints=datapoints
    if mode == :linear
        params.radiusrange = range(start, radius, length=datapoints)
        params.radius = collect(params.radiusrange)
    elseif mode == :log10
        runit = unit(typeof(radius))
        istart = uconvert(u"kpc",start)
        iradius = uconvert(u"kpc",radius)
        params.radiusrange = (range(log10(ustrip(istart)), log10(ustrip(iradius)), length=datapoints))
        params.radius = (10 .^collect(params.radiusrange))u"kpc"
        params.radius = uconvert.(runit, params.radius)    
        istart = uconvert(runit,start)
        iradius = uconvert(runit,radius)
    end

    println("range: R=", params.radiusrange, " ($(params.mode))")
    if mode == :log10
    println("correspond to: R=$(istart) - $(iradius) (linear)")
    end
    
    println("datapoints: ", params.datapoints)
    return params
end