# dark matter
# =====================
function m_profile_NFW(ρc, Rs , r)
    ρc = (ρc)u"Msun/pc^3"
    r = uconvert(u"pc", r)
    Rs = uconvert(u"pc", Rs)

    ρc = ustrip(ρc)
    r = ustrip(r)
    Rs = ustrip(Rs)
    integral = QuadGK.quadgk.(r -> 4 * pi * r^2 * ρ_profile_NFW(ρc, Rs, r), 0. , r, rtol=1e-6,order=7)
    return integral[1]
end

function ρ_profile_NFW(ρc, Rs, r) 
    ρc / (  (r/Rs) * (1 + (r/Rs) )^2  )
end

function fMtot_NFW(ρc, Mtot, Rs, r) 
    return ustrip(Mtot) - m_profile_NFW(ρc, Rs, r)
end  

function v_profile_NFW(ρc, Rs, r)
    sqrt( G * m_profile_NFW(ρc, Rs , r)u"Msun" / r )
end

function Σ_profile_NFW(ρc, Rs, zheight, R)
    ρc = (ρc)u"Msun/pc^3"
    R = uconvert(u"pc", R)
    Rs = uconvert(u"pc", Rs)
    zheight = uconvert(u"pc", zheight)
    
    ρc = ustrip(ρc)
    R = ustrip(R)
    Rs = ustrip(Rs)
    zheight = ustrip(zheight)
    
    integral = QuadGK.quadgk.(z -> ρc / ( sqrt(R^2 + z^2)/Rs * (1. + (sqrt(R^2 + z^2)/Rs))^2 ), -zheight , zheight, rtol=1e-12,order=7)
    return integral[1]u"Msun/pc^2"
end

