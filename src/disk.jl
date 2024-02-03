# constant σ or cs/T
# =====================
function zh_profile(Σc, h, R, σ)
    σ^2 / (Σc * pi * G)  * exp(R/h)
end

function ρ_midplane(ρc, h, R)
    ρc * exp(-R/h)^2
end

# pressure gradient correction
function v_grad(h, R, σ, γ)
    σ * sqrt( 2. * R / (h * γ))
end


function zhc_find(Σc, σ)
    σ^2/(pi * G * Σc)
end


# constant scale height
# =====================
function σc_find(Σc, zh)
    sqrt(zh * pi * G * Σc)
end

function σ_profile(σc, h, R)
    σc * sqrt(exp(-R/h))
end

function ρ_midplane_zh(ρc, h, R)
    ρc * exp(-R/h)
end

# pressure gradient correction
function v_grad_zh(h, R, σ, γ)
    σ * sqrt( R / (h * γ))
end


# exponential disk
# ================
function Σ_profile(Σc, h, R)
    Σc * exp(-R/h)
end

function fMtot(Σc, gset, Mtot, h, R) 
    Σc = (Σc)u"Msun/pc^2"
    R = uconvert(u"pc", R)
    return Mtot - m_profile(Σc, h, R)
end  

function m_profile(Σc, h, R)
    2. * pi * Σc *  h * (h - exp(-R/h) * (h + R) ) # integrated exponential Σ-profile
end

function ρc_find(Σc, σ)
    Σc^2 * pi * G / (2. * σ^2)
end


function y(R, h)
    uconvert(NoUnits, R / (2. * h))
end

function v_profile(Σc, h, R)
    sqrt( 4. * pi * G * Σc * h * y(R,h)^2 * ( besseli(0, y(R,h)) * besselk(0, y(R,h)) -  besseli(1, y(R,h)) * besselk(1, y(R,h)) ) )
end



function T_convert(cs, μ, γ)
    cs^2 * μ * mH / (γ * kB)
end

function cs_convert(T, μ, γ)
    sqrt(γ * kB * T / (μ * mH))
end


