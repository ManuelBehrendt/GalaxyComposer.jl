function overview(params; verbose=true)
    param = Dict()
    profile = Dict()
    
    
    if verbose
        println("Overview (rounded values on screen):")
        println("range: R=", params.radiusrange, " ($(params.mode))")
    end

    if params.mode == :log10
        istart = (10^params.radiusrange[1])u"kpc";      istart_p = runit(istart, 2)
        iradius = (10^params.radiusrange[end])u"kpc";   iradius_p = runit(iradius, 2)
        if verbose println("correspond to: R=$(istart_p) - $(iradius_p) (linear)") end

        Rend = iradius
    else
        Rend = params.radiusrange[end]
    end
    
    if verbose println("datapoints: ", params.datapoints) end
    if length(params.components) ==  0
        if verbose println("components: none set") end
    else
        if verbose println("components: ", keys(params.components) ) end
    end

    
    for ii in keys(params.components)
        comp_param = Dict()
        comp_prof = Dict()
        if verbose println() end
        i = params.components[ii]
        
        
        if i.component == :gasdisk || i.component == :stellardisk
            if verbose printstyled(ii, "\n"; bold = :true) end
            prof = i.profiles
            st = i.structure

            comp_param[:Radius] = (0.0u"kpc", i.R0, Rend)
            if verbose
                R0_p = runit(comp_param[:Radius][2], 3);    Rmax_p = runit(comp_param[:Radius][3], 3)
                println("R0=",R0_p, ", Rmax=",Rmax_p) 
            end

            comp_param[:Mass] = (0.0u"Msun", i.Mtot, prof.Mcum[end])
            if verbose
                MR0_p = runit(comp_param[:Mass][2],3) ; MRmax_p = runit(comp_param[:Mass][3],3)
                println("M(R0)=",MR0_p, ", M(Rmax)=", MRmax_p) 
            end
            if length(prof.V) > 2
                comp_param[:Vrot] = (0.0u"km/s", uconvert(u"km/s", v_profile(st.Σc, st.h, i.R0) ), prof.V[end]) 
                if verbose
                    VR0_p = runit(comp_param[:Vrot][2], 3) ; VRmax_p = runit(comp_param[:Vrot][3], 3)
                    println("V(R0)=",VR0_p , ", V(Rmax)=", VRmax_p ) 
                end
            end
            if length(prof.Vgrad) > 2
                if st.zh != 0.0u"kpc"
                    σR0 = σ_profile(st.σc, st.h, i.R0)
                    comp_param[:Vgrad] = (0.0u"km/s", uconvert(u"km/s", v_grad_zh(st.h, i.R0, σR0, i.γ) ), prof.Vgrad[end])
                else
                    comp_param[:Vgrad] = (0.0u"km/s", uconvert(u"km/s", v_grad(st.h, i.R0, st.σ, i.γ) ), prof.Vgrad[end])
                    
                end
                if verbose
                    VgradR0_p = runit(comp_param[:Vgrad][2],3); VgradRmax_p = runit(comp_param[:Vgrad][3],3)
                    println("Vgrad(R0)=",VgradR0_p , ", Vgrad(Rmax)=", VgradRmax_p )
                end
            end
                
            if i.component == :gasdisk
                comp_param[:μ] = i.μ
                comp_param[:γ] = i.γ
                if verbose println("μ= $(i.μ), γ= $(round(i.γ,sigdigits=3)) (needed for ρ, zh, Vgrad, cs, T)") end
            end

            if verbose
                println()
                println("Details:")
            end
            comp_param[:h] = st.h
            if verbose println("h=", runit(st.h,3), " | scale length") end
            if st.zh != 0.0u"kpc" # constant scale height
                comp_param[:zh] = st.zh
                if verbose printstyled("zh=", runit(comp_param[:zh],3), " | constant scale height", "\n"; color = :red) end
                #println("zh=", comp_param[:zh], " | constant scale height")
                σR0 = σ_profile(st.σc, st.h, i.R0)
                σRmax = σ_profile(st.σc, st.h, Rend)
                comp_param[:σ] = (st.σc, uconvert(u"km/s", σR0), uconvert(u"km/s", σRmax))
                if verbose 
                    σc_p = runit(comp_param[:σ][1], 3) ; σR0_p = runit(comp_param[:σ][2], 3) ; σRmax_p = runit(comp_param[:σ][3], 3)
                    println("σc=",σc_p ,", σ(R0)=",σR0_p , ", σ(Rmax)=",σRmax_p ," | velocity dispersion profile") 
                end

                comp_param[:T] = (st.Tc, uconvert(u"K", T_convert(σR0, i.μ, i.γ)), uconvert(u"K", T_convert(σRmax, i.μ, i.γ)))
                if verbose
                    Tc_p = runit(comp_param[:T][1],3); TR0_p = runit(comp_param[:T][2],3); TRmax_p = runit(comp_param[:T][3],3)
                    println("Tc=",Tc_p ,", T(R0)=",TR0_p , ", T(Rmax)=",TRmax_p ," | temperature profile") 
                end
            else
                comp_param[:zh] = (st.zhc, uconvert(u"pc", zh_profile(st.Σc, st.h, i.R0, st.σ)), uconvert(u"pc", zh_profile(st.Σc, st.h, Rend, st.σ)))
                if verbose
                    zhc_p = runit(comp_param[:zh][1],3); zhR0_p = runit(comp_param[:zh][2],3); zhRmax_p = runit(comp_param[:zh][3],3)
                    println("zhc=",zhc_p ,", zh(R0)=",zhR0_p , ", zh(Rmax)=",zhRmax_p ," | scale height profile") 
                end
                
                comp_param[:σ] = st.σ
                
                if verbose printstyled("σ=", runit(comp_param[:σ], 3), " | constant velocity dispersion", "\n"; color = :red) end
                #println("σ=",comp_param[:σ], " | constant velocity dispersion")

                comp_param[:T] = st.T
                if verbose printstyled("T=", runit(comp_param[:T], 3),  " | constant temperature", "\n"; color = :red) end
                #println("T=",comp_param[:T],  " | constant temperature")
            end

            comp_param[:Σ] = (st.Σc, uconvert(u"Msun/pc^2", Σ_profile(st.Σc, st.h, i.R0)) ,  Σ_profile(st.Σc, st.h, Rend))
            if verbose
                Σc_p = runit(comp_param[:Σ][1],3); ΣR0_p = runit(comp_param[:Σ][2],3); ΣRmax_p = runit(comp_param[:Σ][3],3)
                println("Σc=",Σc_p , ", Σ(R0)=",ΣR0_p , ", Σ(Rmax)=",ΣRmax_p ," | surface density profile") 
            end
            if st.ρc != 0.0u"Msun/pc^3"
                if st.zh != 0.0u"kpc"
                    comp_param[:ρ] = (st.ρc, uconvert(u"Msun/pc^3", ρ_midplane_zh(st.ρc, st.h, i.R0)),  uconvert(u"Msun/pc^3", ρ_midplane_zh(st.ρc, st.h, Rend)))
                    
                else
                    comp_param[:ρ] = (st.ρc, uconvert(u"Msun/pc^3", ρ_midplane(st.ρc, st.h, i.R0)),  uconvert(u"Msun/pc^3", ρ_midplane(st.ρc, st.h, Rend)))
                end
                
                if verbose
                    ρc_p = runit(comp_param[:ρ][1],3); ρ0R0_p = runit(comp_param[:ρ][2],3); ρ0Rmax_p = runit(comp_param[:ρ][3] ,3)
                    println("ρc=",ρc_p , ", ρ0(R0)=",ρ0R0_p , ", ρ0(Rmax)=",ρ0Rmax_p  ," | midplane density profile")
                end
            end


            comp_prof[:R]     = ustrip(params.radius)
            comp_prof[:Mcum]  = ustrip(prof.Mcum)
            comp_prof[:V]     = ustrip(prof.V)
            comp_prof[:Σ]     = ustrip(prof.Σ)
            if length(prof.Vgrad) > 2
                comp_prof[:Vgrad] = ustrip(prof.Vgrad)
            end
            comp_prof[:zh]    = ustrip(prof.zh)
            comp_prof[:ρ0]    = ustrip(prof.ρ0)
            comp_prof[:σ]     = ustrip(prof.σ)
            comp_prof[:T]     = ustrip(prof.T)

        elseif i.component == :NFW
            if verbose printstyled(ii, "\n"; bold = :true) end
            prof = i.profiles
            st = i.structure

            comp_param[:Radius] = (0.0u"kpc", i.R0, i.Rvir)
            if verbose
                R0_p = runit(comp_param[:Radius][2],3); Rvir_p = runit(comp_param[:Radius][3],3)
                println("R0=",R0_p , ", Rvir=",Rvir_p) 
            end

            comp_param[:Mass] = (0.0u"Msun", i.Mtot, prof.Mcum[end], i.Mvir )
            if verbose 
                MR0_p = runit(comp_param[:Mass][2],3); MRmax_p = runit(comp_param[:Mass][3],3); Mvir_p = runit(comp_param[:Mass][4],3)
                println("M(R0)=",MR0_p , ", M(Rmax)=",MRmax_p , ", Mvir=",Mvir_p  ) 
            end

            if verbose 
                println()
                println("Details:")
            end
            comp_param[:Rs] = st.Rs
            if verbose 
                println("Rs=", runit(st.Rs,3), " | scale radius")
                println("c=", st.c, " | concentration")
            end

            comp_param[:Σ] = (prof.Σ[1], uconvert(u"Msun/pc^2", Σ_profile_NFW(st.ρc, st.Rs, i.zheight, i.R0)) ,  prof.Σ[end])
            if verbose
                ΣRmin_p = runit(comp_param[:Σ][1],3); ΣR0_p = runit(comp_param[:Σ][2],3); ΣRvir_p = runit(comp_param[:Σ][3],3)
                println("Σ(Rmin)=",ΣRmin_p , ", Σ(R0)=",ΣR0_p , ", Σ(Rvir)=",ΣRvir_p ," | surface density profile") 
            end
              
            comp_param[:Vrot] = (0.0u"km/s", uconvert(u"km/s", v_profile_NFW(st.ρc, st.Rs, i.R0) ), prof.V[end]) 
            if verbose
                VR0_p = runit(comp_param[:Vrot][2],3); VRvir_p = runit(comp_param[:Vrot][3],3)
                println("V(R0)=",VR0_p , ", V(Rvir)=",VRvir_p  ) 
            end
            

            comp_prof[:R]     = ustrip(params.radius)
            comp_prof[:Mcum]  = ustrip(prof.Mcum)
            comp_prof[:V]     = ustrip(prof.V)
            comp_prof[:Σ]     = ustrip(prof.Σ)
            comp_prof[:ρ]     = ustrip(prof.ρ)
            
        elseif i.component == :bulge
            if verbose printstyled(ii, "\n"; bold = :true) end
            prof = i.profiles
            st = i.structure

            comp_param[:Radius] = (0.0u"kpc", i.R0, Rend)
            if verbose
                R0_p = runit(comp_param[:Radius][2], 3);    Rmax_p = runit(comp_param[:Radius][3], 3)
                println("R0=",R0_p, ", Rmax=",Rmax_p) 
            end


            comp_param[:Mass] = (0.0u"Msun", i.Mtot, prof.Mcum[end])
            if verbose
                MR0_p = runit(comp_param[:Mass][2],3) ; MRmax_p = runit(comp_param[:Mass][3],3)
                println("M(R0)=",MR0_p, ", M(Rmax)=", MRmax_p) 
            end
            if length(prof.V) > 2
                comp_param[:Vrot] = (0.0u"km/s", uconvert(u"km/s", v_profile_bulge(st.ρc, st.a, i.R0) ), prof.V[end]) 
                if verbose
                    VR0_p = runit(comp_param[:Vrot][2], 3) ; VRmax_p = runit(comp_param[:Vrot][3], 3)
                    println("V(R0)=",VR0_p , ", V(Rmax)=", VRmax_p ) 
                end
            end

            comp_param[:Σ] = (st.Σc, uconvert(u"Msun/pc^2", Σ_profile_bulge(st.ρc, st.a, i.R0)) ,  Σ_profile_bulge(st.ρc, st.a, Rend))
            if verbose
                Σc_p = runit(comp_param[:Σ][1],3); ΣR0_p = runit(comp_param[:Σ][2],3); ΣRmax_p = runit(comp_param[:Σ][3],3)
                println("Σc=",Σc_p , ", Σ(R0)=",ΣR0_p , ", Σ(Rmax)=",ΣRmax_p ," | surface density profile") 
            end


            comp_param[:ρ] = (st.ρc, uconvert(u"Msun/pc^3", ρ_profile_bulge(st.ρc, st.a, i.R0)) ,  ρ_profile_bulge(st.ρc, st.a, Rend))
            if verbose
                ρc_p = runit(comp_param[:ρ][1],3); ρ0R0_p = runit(comp_param[:ρ][2],3); ρ0Rmax_p = runit(comp_param[:ρ][3] ,3)
                println("ρc=",ρc_p , ", ρ0(R0)=",ρ0R0_p , ", ρ0(Rmax)=",ρ0Rmax_p  ," | midplane density profile")
            end

            comp_prof[:R]     = ustrip(params.radius)
            comp_prof[:Mcum]  = ustrip(prof.Mcum)
            comp_prof[:V]     = ustrip(prof.V)
            comp_prof[:Σ]     = ustrip(prof.Σ)
            comp_prof[:ρ]    = ustrip(prof.ρ)

        end
        param[ii] = comp_param
        profile[ii] = comp_prof
    end
    return param, profile
end



function runit(value, digits)
    round(typeof(value), value, sigdigits=digits)
end