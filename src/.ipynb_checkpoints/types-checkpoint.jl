@with_kw mutable struct SetupParams
    radius::Vector{Float64}
    components::Dict{Any, Any}
end

@with_kw mutable struct DiskStructure
    h::Real=0.    # exponential disc scale-length
    Σc::Real=0.   # central surface density
    ρc::Real=0.   # central density
    
    # used for constant scale-height
    zh::Real=0. # constant scale-height
    σc::Real=0. # central velocity dispersion/sound speed  
    Tc::Real=0. # central temperature

    # used for constant "temperature"
    T::Real=0. # constant temperature
    σ::Real=0. # constant velocity dispersion/sound speed  
end

@with_kw mutable struct InitDisk
    component::Symbol=:gas # :gas or :stars
    Mtot::Real=0. # total mass within radius range
    structure::DiskStructure=DiskStructure()
    μ::Real=1.4   # molecular weight
    γ::Real=5/3   # adiabatic index
end


