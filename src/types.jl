
G  = 6.67384e-8u"cm^3/g /s^2"          
kB = 1.380649e-16u"cm^2 * g /s ^2 / K"
mH = 1.0078250322 * 1.66053906660e-27u"kg"

@with_kw mutable struct SetupParams
    radius=collect(1e-3:0.1:16.)u"kpc"
    radiusrange=(1e-3:0.1:16.)u"kpc"
    datapoints::Int=100
    mode::Symbol=:linear
    components::Dict{Any, Any}=Dict()
end


# ======= disk =====================
@with_kw mutable struct DiskStructure
    h=0.0u"kpc"    # exponential disc scale length
    Σc=0.0u"Msun/pc^2"   # central surface density
    ρc=0.0u"Msun/pc^3"   # central density
    
    # used for constant scale height
    zh=0.0u"pc" # constant scale height
    σc=0.0u"km/s" # central velocity dispersion/sound speed  
    Tc=0.0u"K" # central temperature

    # used for constant "temperature"
    zhc=0.0u"pc" # central scale height
    T=0.0u"K" # constant temperature
    σ=0.0u"km/s" # constant velocity dispersion/sound speed  
end

@with_kw mutable struct DiskProfiles
    Mcum=[0.,1.]u"Msun"
    V=[0.,1.]u"km/s"
    Σ=[0.,1.]u"Msun/pc^2"
    Vgrad=[0.,1.]u"km/s"
    zh=[0.,1.]u"pc"
    ρ0=[0.,1.]u"Msun/pc^3"
    σ=[0.,1.]u"km/s"
    T=[0.,1.]u"K"
end

@with_kw mutable struct InitDisk
    component::Symbol=:gas # :gasdisk or :stellardisk
    Mtot=0.0u"Msun" # total mass within R0 range
    R0=0.0u"kpc" # disk radius used for calculations, e.g. Mtot
    structure::DiskStructure=DiskStructure()
    profiles::DiskProfiles=DiskProfiles()
    μ::Real=1.4   # molecular weight
    γ::Real=5/3   # adiabatic index
    potential=false # if used as potential, e.g. no Vgrad is calculated
end



# ======= dark matter =====================
@with_kw mutable struct NFWStructure
    ρc=0.0u"Msun/pc^3" # central density
    Rs=0.0u"kpc" # scale radius
    c::Real=0.  # concentration parameter (dominates over scale radius)
end

@with_kw mutable struct NFWProfiles
    V=[0.,1.]u"km/s"
    ρ=[0.,1.]u"Msun/pc^3"
    Mcum=[0.,1.]u"Msun"
    Σ=[0.,1.]u"Msun/pc^2"
end

@with_kw mutable struct InitNFW
    component::Symbol=:NFW
    Mvir=0.0u"Msun" # virial mass
    Rvir=0.0u"kpc" # virial radius (dominates over scale radius)
    Mtot=0.0u"Msun" # total mass within R0 range
    R0=0.0u"kpc" # radius used for calculations, e.g. Mtot
    structure::NFWStructure=NFWStructure()
    profiles::NFWProfiles=NFWProfiles()
    zheight=0.0u"pc" # integration height for surface density
    potential=true # if used as potential
end

# ======= stellar bulge =====================
@with_kw mutable struct Bulgetructure
    Σc=0.0u"Msun/pc^2"   # central surface density
    ρc=0.0u"Msun/pc^3" # central density
    a=0.0u"pc" # scale radius
end

@with_kw mutable struct BulgeProfiles
    V=[0.,1.]u"km/s"
    ρ=[0.,1.]u"Msun/pc^3"
    Mcum=[0.,1.]u"Msun"
    Σ=[0.,1.]u"Msun/pc^2"
end

@with_kw mutable struct InitBulge
    component::Symbol=:bulge
    Mtot=0.0u"Msun" # total mass within R0 range
    R0=0.0u"kpc" # radius used for calculations, e.g. Mtot
    structure::Bulgetructure=Bulgetructure()
    profiles::BulgeProfiles=BulgeProfiles()
    potential=true # if used as potential
end