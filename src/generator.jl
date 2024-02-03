

function generate(gset::SetupParams, component::Symbol, number::Int; kwargs...)
 
    ks = keys(kwargs)

    verbose = false
    if :verbose in ks verbose = kwargs[:verbose] end
    if verbose print("\n generate: ", Symbol(String(component) * "$number")) end

    if component == :gasdisk || component == :stellardisk
        if :Mtot in ks && :h in ks
    
            disk = InitDisk()
            disk.component = component
            if component == :stellardisk disk.potential = true end
            if :potential in ks disk.potential = kwargs[:potential] end
            if verbose && disk.potential == true print("\t for use as potential") end
            
            # use familiar units
            Mtot = uconvert(u"Msun", kwargs[:Mtot])
            h = uconvert(u"pc", kwargs[:h])
            if :R0 in ks
                R0 = kwargs[:R0]
            else
                R0 = gset.radius[end]
            end
            
            # find Σc for given Mtot, h at R0
            x = find_zero(x ->fMtot(x, gset, Mtot, h, R0), 1e2)
            Σc = (x)u"Msun/pc^2"

            # remember parameter
            disk.Mtot = Mtot
            disk.structure.h = uconvert(u"kpc", h)
            disk.structure.Σc = Σc
            disk.R0 = uconvert(u"kpc",R0)

            # get profiles
            disk.profiles.Mcum = @. uconvert(u"Msun", m_profile(Σc, h, gset.radius))
            disk.profiles.V =  @. uconvert(u"km/s", v_profile(Σc, h, gset.radius))
            disk.profiles.Σ =  @. uconvert(u"Msun/pc^2", Σ_profile(Σc, h, gset.radius))

            
            if :σ in ks  # assuming constant σ -> flaring disk
                σ = kwargs[:σ]
                gen_const_σ(gset, disk, Σc, h, σ)
                
            elseif :T in ks # assuming constant T -> flaring disk
                T = kwargs[:T]
                gen_const_T(gset, disk, Σc, h, T)
  

            elseif :zh in ks # assuming constang scale height -> decreasing σ or T
                zh = kwargs[:zh]
                gen_const_zh(gset, disk, Σc, h, zh)

            end
            gset.components[Symbol(String(component) * "$number")] = disk

        end # disk

    elseif component == :NFW
        if :Mtot in ks 
            dm = InitNFW()
            dm.component = component
            #if :potential in ks dm.potential = kwargs[:potential] end
            if verbose && dm.potential == true print("\t for use as potential") end
            
            Mtot = uconvert(u"Msun", kwargs[:Mtot])
            
            if :c in ks && :Rvir in ks
                c = kwargs[:c]
                Rvir = kwargs[:Rvir]
                Rs = Rvir / c
            elseif :c in ks && :Rs in ks
                c = kwargs[:c]
                Rs = kwargs[:Rs]
                Rvir = Rs * c
            elseif :Rvir in ks && :Rs in ks
                Rvir = kwargs[:Rvir]
                Rs = kwargs[:Rs]
                c = Rvir / Rs
            end


            dm.Rvir = uconvert(u"kpc", Rvir)
            dm.structure.c = uconvert(NoUnits, c)
            dm.structure.Rs = uconvert(u"kpc", Rs)

            if :R0 in ks
                R0 = kwargs[:R0]
            else
                R0 = dm.Rvir
            end


            # find ρc for given Mtot, Rs at R0
            x = find_zero(x ->fMtot_NFW(x, Mtot, Rs, R0), 0.1)
            ρc = (x)u"Msun/pc^3"
            Mvir = m_profile_NFW(ρc, Rs, Rvir)u"Msun"
            dm.Mvir = Mvir
            
            # remember parameter
            dm.Mtot = Mtot
            dm.structure.ρc = ρc
            dm.R0 = uconvert(u"kpc", R0)

            # vertical integration range vor surface density
            if :zheight in ks
                zheight = kwargs[:zheight]
            else
                zheight = Rvir
            end
            dm.zheight = uconvert(u"pc", zheight)
            
            # get profiles
            dm.profiles.Mcum = @. m_profile_NFW(ρc, Rs, gset.radius)u"Msun"
            dm.profiles.V =  @. uconvert(u"km/s", v_profile_NFW(ρc, Rs, gset.radius))
            dm.profiles.ρ =  @. uconvert(u"Msun/pc^3", ρ_profile_NFW(ρc, Rs, gset.radius))
            dm.profiles.Σ =  @. Σ_profile_NFW(ρc, Rs, zheight, gset.radius)
        end # NFW


        gset.components[Symbol(String(component) * "$number")] = dm

    elseif component == :bulge 
        if :Mtot in ks && :a in ks
            bulge = InitBulge()
            bulge.component = component
            #if :potential in ks bulge.potential = kwargs[:potential] end
            if verbose && bulge.potential == true print("\t for use as potential") end 
            
            Mtot = uconvert(u"Msun", kwargs[:Mtot])

            a = kwargs[:a]
            bulge.structure.a = uconvert(u"pc", a)
            
            if :R0 in ks
                R0 = kwargs[:R0]
            else
                R0 = gset.radius[end]
            end

            # find ρc for given Mtot, a at R0
            x = find_zero(x ->fMtot_bulge(x, Mtot, a, R0), 1.)
            ρc = (x)u"Msun/pc^3"
            Mbulge = m_profile_NFW(ρc, a, R0)u"Msun"
            bulge.Mtot = Mbulge
            bulge.structure.ρc = ρc
            bulge.R0 = uconvert(u"kpc", R0)
            bulge.structure.Σc =  Σ_profile_bulge(ρc, a, 0.0u"kpc")

             # get profiles
             bulge.profiles.Mcum = @. m_profile_bulge(ρc, a, gset.radius)u"Msun"
             bulge.profiles.V =  @. uconvert(u"km/s", v_profile_bulge(ρc, a, gset.radius))
             bulge.profiles.ρ =  @. uconvert(u"Msun/pc^3", ρ_profile_bulge(ρc, a, gset.radius))
             bulge.profiles.Σ =  @. Σ_profile_bulge(ρc, a, gset.radius)
        end # Bulge
        gset.components[Symbol(String(component) * "$number")] = bulge
    end
end



function generate(gset::SetupParams, activity::Symbol; verbose = true, components=[:all])
    total = Dict()
    if activity == :total
        if components==[:all]
            icomp = Symbol.(keys(gset.components))
        else
            icomp = deepcopy(components)
        end

    
        Vtot_2 = [0.]  
        Mcum = [0.]
        firsttime = true
        for i in icomp
            data = gset.components[i]
            if verbose println(i) end

            if firsttime
                Vtot_2 = data.profiles.V .^2
                Mcum = data.profiles.Mcum
                firsttime = false
            else
                Vtot_2 .+= data.profiles.V .^2 
                Mcum .+= data.profiles.Mcum
            end

            if data.component == :gasdisk || data.component == :stellardisk
                if length(data.profiles.Vgrad) > 2
                    Vtot_2 .-= data.profiles.Vgrad .^2 
                end
            end
            
        end
        total[:R] = ustrip(gset.radius)
        total[:Vrot] =   ustrip(sqrt.(Vtot_2))
        total[:Ω]     = ustrip(total[:Vrot]) ./ ustrip(gset.radius)
        total[:Mcum] =   ustrip(Mcum)
    end
    return total
end







# generate exponential disk

function gen_ρc(disk, Σc, σ)
    ρc = Σc^2 * pi * G / (2. * σ^2)
    disk.structure.ρc = uconvert(u"Msun/pc^3", ρc)
    return ρc
end


function gen_const_σ(gset, disk, Σc, h, σ)
    disk.profiles.zh =  @. uconvert(u"pc", zh_profile(Σc, h, gset.radius, σ))
    ρc = ρc_find(Σc, σ)
    disk.structure.ρc = uconvert(u"Msun/pc^3", ρc)

    disk.profiles.ρ0 =  @. uconvert(u"Msun/pc^3", ρ_midplane(ρc, h, gset.radius))


    
    disk.structure.σ = uconvert(u"km/s", σ)
    disk.structure.T = uconvert(u"K", T_convert(σ, disk.μ, disk.γ))
    disk.structure.zhc = uconvert(u"pc", zhc_find(Σc, σ))
    disk.profiles.σ = fill(σ, length(gset.radius))u"km/s"
    disk.profiles.T = fill( uconvert(u"K", T_convert(σ, disk.μ, disk.γ)), length(gset.radius))

    if disk.potential == false
        disk.profiles.Vgrad = @. uconvert(u"km/s", v_grad(h, gset.radius, σ, disk.γ))
    end
    return 
end



function gen_const_T(gset, disk, Σc, h, T)

    disk.structure.T = uconvert(u"K", T)
    disk.structure.σ = uconvert(u"km/s", cs_convert(T, disk.μ, disk.γ))
    disk.structure.zhc = uconvert(u"pc", zhc_find(Σc, σ))
    disk.profiles.T = fill( uconvert(u"K", T), length(gset.radius))
    disk.profiles.σ = uconvert.(u"km/s", cs_convert.(disk.profiles.T, disk.μ, disk.γ) )
    σ = disk.profiles.σ
    
    disk.profiles.zh =  @. uconvert(u"pc", zh_profile(Σc, h, gset.radius, σ[1]))
    
    ρc = ρc_find(Σc, σ[1])
    disk.structure.ρc = uconvert(u"Msun/pc^3", ρc)
    
    disk.profiles.ρ0 =  @. uconvert(u"Msun/pc^3", ρ_midplane(ρc, h, gset.radius))
    if disk.potential == false
        disk.profiles.Vgrad = @. uconvert(u"km/s", v_grad(h, gset.radius, σ[1], disk.γ)) 
    end
    return
end


function gen_const_zh(gset, disk, Σc, h, zh)

    disk.structure.zh = uconvert(u"pc",zh)
    disk.profiles.zh =  fill( uconvert(u"pc",zh), length(gset.radius)) 
    
    σc = @. uconvert( u"km/s", σc_find(Σc, zh))
    disk.structure.σc = σc
    disk.structure.Tc = uconvert(u"K", T_convert(σc, disk.μ, disk.γ))
    disk.profiles.σ = @. σ_profile(σc, h, gset.radius)
    disk.profiles.T = @. uconvert(u"K", T_convert(disk.profiles.σ, disk.μ, disk.γ)) 
    
    ρc = ρc_find(Σc, σc)
    disk.structure.ρc = uconvert(u"Msun/pc^3", ρc)
    disk.profiles.ρ0 =  @. uconvert(u"Msun/pc^3", ρ_midplane_zh(ρc, h, gset.radius))
    if disk.potential == false
        disk.profiles.Vgrad = @. uconvert(u"km/s", v_grad_zh(h, gset.radius, disk.profiles.σ, disk.γ)) 
    end
    return
end