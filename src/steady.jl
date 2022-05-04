#=
  Routines for solving the steady-state sausage equations.
=#


"""
   Solve the steady-state sausage model (Luque PSST 2017) to compute a 
   charge distribution λ, head field eh and dimensionless correction to the
   surface current moment, ξ1
"""
function solve_sausage_attach(a, v, eb; ne0=1e20, τatt=40e-9, L=0.03, Δzmin=1e-5,
                              nlevs=5, Δzeaxis=5e-7)

    zf = refined_mesh(L, Δzmin, nlevs)
    zc = @. 0.5 * (zf[2:end] + zf[1:end - 1])
    Δz = diff(zf)
    
    # Last edge removed.
    zf1 = zf[1:end-1]
    
    MC = zeros(length(zf) - 1, length(zc))
    MS = zeros(length(zf) - 1, length(zc))
    MD = zeros(length(zf) - 1, length(zc))
    
    ξ = 4.0

    ϵ = 1e-5

    latt = v * τatt
    μe = 0.0372193
    σ = @. co.e * ne0 * μe * exp(zf1 / latt)
    
    δ = ξ * v * co.ϵ0 / σ[end]

    A = @. σ / (4π * co.ϵ0)
    @info "The matrix size is" size(MC)
    
    Threads.@threads for j in reverse(eachindex(zc))
        z1 = zc[j]
        for (i, z) in enumerate(zf1)
            Δz = zf[j + 1] - zf[j]
            MC[i, j] = A[i] * GC_GK(a, z, zf[j], zf[j + 1])
            MS[i, j] = A[i] * GS_GK(a, ϵ, δ, z, zf[j], zf[j + 1])

            if i == j
                MD[i, j] = -0.5 * v
                # Neumann b.c. at -L
                if i == 1
                    MD[i, j] -= 0.5 * v
                end
            end

            if i == j + 1
                MD[i, j] = -0.5 * v
            end                
        end
    end    
    M = MC + MS + MD

    @info "Matrix assembly completed"
    
    I0 = @. σ * π * eb * R(a, zf1)^2 

    λ = -M \ I0
    
    eh = eb
    el = eb

    #Evaluation of the full axial field
    zaxis = -L/2:Δzeaxis:L/2
    ME = zeros(length(zaxis), length(zc))

    # A small distance to prevent computeing stuff just at the discontinuity
    δz = 1e-9
    Threads.@threads for j in eachindex(zc)
        for (i, z) in enumerate(zaxis)
            #ME[i, j] = Δz / (4π * co.ϵ0) * GR(a, z, z1)
            # Slight shift to prevent Inf inside the integral.
            zs = z == 0 ? δz : z
            ME[i, j] = 1 / (4π * co.ϵ0) * GR_GK(a, zs, zf[j], zf[j + 1])
        end
        eh += (1 / (4π * co.ϵ0)) * λ[j] * GR_GK(a,  δz, zf[j], zf[j + 1])
        el += (1 / (4π * co.ϵ0)) * λ[j] * GR_GK(a, -δz, zf[j], zf[j + 1])        
    end
    
    eaxis = eb .+ ME * λ

    (;λ, eh, el, zaxis, eaxis, σ, zc, zf)
end


function refined_mesh(L, Δzmin, nlevs)
    s = L / (2^nlevs)
    n = round(Int, s / Δzmin)
    
    zf = range(-s, 0, length=n + 1)

    for l in 1:nlevs
        s *= 2
        zf = [range(-s, -s/2, length=n + 1); zf[2:end]]
    end
    
    return zf
end
