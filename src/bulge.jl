
function ρ_profile_bulge(ρc, a, R) 
    ρc * exp(- R / a)
end


function m_profile_bulge(ρc, a, R)
    ρc = (ρc)u"Msun/pc^3"
    R = uconvert(u"pc", R)
    a = uconvert(u"pc", a)

    ρc = ustrip(ρc)
    R = ustrip(R)
    a = ustrip(a)
    integral = QuadGK.quadgk.(R -> 4 .* pi .* R^2 .* ρ_profile_bulge(ρc, a, R), 0. , R, rtol=1e-6,order=7)
    return integral[1]
end

function fMtot_bulge(ρc, Mtot, a, R) 
    return ustrip(Mtot) - m_profile_bulge(ρc, a, R)
end  

function v_profile_bulge(ρc, a, R)
    sqrt( G * m_profile_bulge(ρc, a , R)u"Msun" / R )
end

function Σ_profile_bulge(ρc, a, R)
    ρc = (ρc)u"Msun/pc^3"
    R = uconvert(u"pc", R)
    a = uconvert(u"pc", a)
    
    ρc = ustrip(ρc)
    R = ustrip(R)
    a = ustrip(a)
    
    integral = QuadGK.quadgk.(z -> ρc * exp(-sqrt(R^2 + z^2)/a), -100. , 100., rtol=1e-12,order=7)
    return integral[1]u"Msun/pc^2"
end
