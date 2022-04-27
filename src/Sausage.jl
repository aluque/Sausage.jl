module Sausage

using QuadGK
using Constants: co
using LinearAlgebra
using Roots
using SpecialFunctions

include("naidis.jl")

R(a, z) = a * sqrt(1 - exp(z / a))
xf(a, z1, ϕ1) = R(a, z1) * cos(ϕ1)
ρf(a, z, z1, ϕ1) = sqrt(R(a, z1)^2 * sin(ϕ1)^2 + (z - z1)^2)

function GC1(a, z, z1, ϕ1, r)
    ρ = ρf(a, z, z1, ϕ1)
    x = xf(a, z1, ϕ1)

    (x * (r - x) - ρ^2) / (ρ^2 * sqrt(ρ^2 + (r - x)^2))
end
    
function GC(a, z, z1)
    (z ≈ z1) && return 0.0 
    (z - z1) * quadgk(ϕ1 -> (GC1(a, z, z1, ϕ1, R(a, z)) - GC1(a, z, z1, ϕ1, 0)),
                      0, 2π)[1]
end


function U1(a, sϵ, z, z1, ϕ1)
    ρ = ρf(a, z, z1, ϕ1)
    x = xf(a, z1, ϕ1)

    1 / (2π * ((R(a, z) + sϵ - x)^2 + ρ^2)^(3/2))
end

function U(a, sg, ϵ, z, z1)
    (z - z1) * quadgk(ϕ1 -> U1(a, sg * ϵ, z, z1, ϕ1), 0, 2π)[1]
end

function GS(a, ϵ, δ, z, z1)
    (z ≈ z1) && return 0.0 

    R⁻ = R(a, z)
    R⁺ = R(a, z - δ)

    (π / 6) * (R⁺ - R⁻) * (U(a, +1, ϵ, z, z1) * (R⁺ - R⁻) +
                           U(a, -1, ϵ, z, z1) * (R⁺ + 3R⁻))
end

function GR(a, z, z1)
    (z - z1) / ((R(a, z1)^2 + (z - z1)^2)^(3/2))
end

function solve_velocity_a(eb, a)
    vmin, vmax = 1e5, 1e8
    println()
    @show eb
    
    v0 = find_zero((vmin, vmax), Roots.A42()) do v
        λ, eh = solve_sausage(a, v, eb)
        f = v - velocity(+1, eh, a)
    end
    @show v0
    v0
    #(;λ, eh)
end

function solve_velocity_eh(eb, eh)
    println()
    @show eb

    amin, amax = 1e-4, 5e-3
    
    a = find_zero((amin, amax), Roots.A42()) do a
        v = velocity(+1, eh, a)
        @show a v

        λ, eh1 = solve_sausage(a, v, eb)
        @show eh1 eh
        
        f = eh1 - eh
    end

    velocity(+1, eh, a / 2.5)
end

function solve_velocity_briels(eb)
    vmin, vmax = 1e5, 5e7
    println()
    @show eb
    c = 2e12 / 4

    v = find_zero((vmin, vmax), Roots.A42()) do v1
        a = sqrt(v1 / c)
        λ, eh = solve_sausage(a, v1, eb)
        @show a v1 eh
        f = v1 - velocity(+1, eh, a)
    end

    a = sqrt(v / c)
    (;v, a)
    
    #(;λ, eh)
end

function solve_radius_naidis(eb)
    println()
    @show eb

    amin, amax = 1e-4, 15e-3
    va = vallen(eb)

    a = find_zero((amin, amax), Roots.A42()) do a
        λ, eh = solve_sausage(a, va, eb)

        vn = velocity(+1, eh, a)
        @show a (vn - va)
        f = vn - va
    end

    a
end


"""
    Finds a set of eb, eh, v, a, ne that satisfy all the following constraints
    - Naidis
    - Allen & Mikropoulos
    - Lagarkov (ionization integral)
    - Sausage steady state.

    To avoid solving more than one unknown we start from the head field although
    we typically want everything as a function of the background field.
"""    
function solve_full(eh)
    # We can directly obtain ne
    ne = ne_integral(eh)

    # we solve for a in this range
    amin, amax = 1e-4, 15e-3

    a = find_zero((amin, amax), Roots.A42(), xrtol=1e-4) do a
        @show a
        
        v, eb, eh1 = solve_inner(eh, a, ne)
        
        # the resulting eh has to equal the given one
        eh1 - eh
    end

    v, eb, eh1 = solve_inner(eh, a, ne)

    (;eb, v, a, ne)
end


function solve_inner(eh, a, ne)
    # at this point we have eh, a.  With naidis we obtain v
    v = velocity(+1, eh, a)
    
    # inverting Allen & Mikropoulos we get eb
    eb = inv_vallen(v)
    
    # now we solve for the steady-state in the sausage model
    _, eh1, _ = solve_sausage(a, v, eb, ne=ne)

    (;v, eb, eh1)
end


"""
   Solve the system that includes
    - Naidis
    - Lagarkov (ionization integral)
    - Sausage steady state.
    - Briels
   We start with the velocity to avoid solveing 2d equations (but we solve
   2 consecutive 1d equations.
"""
function solve_noallen(v)
    # Briels
    c = 2e12 / 4
    a = sqrt(v / c)

    # solve for Naidis
    ehmin, ehmax = (5e6, 3e7)
    eh = find_zero((ehmin, ehmax), Roots.A42(), xrtol=1e-4) do eh
        v - velocity(+1, eh, a)
    end

    # ionization
    ne = ne_integral(eh)

    # solve for Sausage
    ebmin, ebmax = (1e5, 5e6)
    eb = find_zero((ebmin, ebmax), Roots.A42(), xrtol=1e-4) do eb
        _, eh1, _ = solve_sausage(a, v, eb)
        eh1 - eh
    end
    (;eb, eh, a, ne)

end


"""
   Solve the steady-state sausage model (Luque PSST 2017) to compute a 
   charge distribution λ, head field eh and dimensionless correction to the
   surface current moment, ξ1
"""
function solve_sausage(a, v, eb; ne=1e20)
    L = 1.0
    Δz =20e-5

    ξ = 4.0

    ϵ = 1e-5

    μe = 0.0372193
    σ = co.e * ne * μe
    
    δ = ξ * v * co.ϵ0 / σ

    zf = -L:Δz:0.0
    zc = @. 0.5 * (zf[2:end] + zf[1:end - 1])
    
    MC = zeros(length(zc), length(zc))
    MS = zeros(length(zc), length(zc))
    
    A = Δz * σ / (4π * co.ϵ0)

    Threads.@threads for j in eachindex(zc)
        z1 = zc[j]
        for (i, z) in enumerate(zc)
            MC[i, j] = A * GC(a, z, z1)
            MS[i, j] = A * GS(a, ϵ, δ, z, z1)
        end
    end    
    M = MC + MS
    
    I0 = @. σ * π * eb * R(a, zc)^2 

    λ = -(M - I * v) \ I0
    
    eh = eb + (Δz / (4π * co.ϵ0)) * sum(λ .* GR.(Ref(a), 0, zc))
    ξ1 = Δz * sum(MS * λ) / (π * a^2 * v * co.ϵ0 * eh)
    
    (;λ, eh, ξ1)
end


"""
   Solve the steady-state sausage model (Luque PSST 2017) to compute a 
   charge distribution λ, head field eh and dimensionless correction to the
   surface current moment, ξ1
"""
function solve_sausage_attach(a, v, eb; ne0=1e20, τatt=70e-9)
    L = 0.05
    Δz =1e-5

    zf = -L:Δz:0.0
    zc = @. 0.5 * (zf[2:end] + zf[1:end - 1])

    # Last edge removed.
    zf1 = zf[1:end-1]
    
    MC = zeros(length(zf) - 1, length(zc))
    MS = zeros(length(zf) - 1, length(zc))
    ME = zeros(length(zf) - 1, length(zc))
    MD = zeros(length(zf) - 1, length(zc))
    
    ξ = 4.0

    ϵ = 1e-5

    latt = v * τatt
    μe = 0.0372193
    σ = @. co.e * ne0 * μe * exp(zf1 / latt)
    
    δ = ξ * v * co.ϵ0 / σ[end]

    A = @. Δz * σ / (4π * co.ϵ0)
    @info "The matrix size is" size(MC)
    
    Threads.@threads for j in eachindex(zc)
        z1 = zc[j]
        for (i, z) in enumerate(zf1)
            MC[i, j] = A[i] * GC(a, z, z1)
            MS[i, j] = A[i] * GS(a, ϵ, δ, z, z1)
            ME[i, j] = Δz / (4π * co.ϵ0) * GR(a, z, z1)

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
    
    eh = eb + (Δz / (4π * co.ϵ0)) * sum(λ .* GR.(Ref(a), 0, zc))
    ech = eb .+ ME * λ

    ξ1 = Δz * sum(MS * λ) / (π * a^2 * v * co.ϵ0 * eh)
    
    (;λ, eh, ξ1, σ, zc, zf)
end


""" 
    Solve the system to satisfy A&M, Briels, Naidis and the Ionization int.
"""
function solve_abni(eb)
    # A&M
    v = vallen(eb)

    # Briels
    c = 2e12 / 4
    a = sqrt(v / c)

    # solve for Naidis
    ehmin, ehmax = (5e6, 3e7)
    eh = find_zero((ehmin, ehmax), Roots.A42(), xrtol=1e-4) do eh
        v - velocity(+1, eh, a)
    end

    # ionization
    ne = ne_integral(eh)

    (;eh, v, a, ne)
end

""" 
    Velocity from Allen & Mikropoulos. 
"""
function vallen(E)
    v0 = 1.25e5
    E0 = 4.91e5
    α = 3.0

    v0 * (E / E0)^α
end


""" 
    Eb from the velocity from Allen & Mikropoulos. 
"""
function inv_vallen(v)
    v0 = 1.25e5
    E0 = 4.91e5
    α = 3.0

    E0 * (v / v0)^(1 / α)
end


""" 
    Ionization integral.
    ∫ exp(-E0/ϵ) dϵ from E1 to E2.

"""
function I0(E0, E1, E2)
    E2 / exp(E0 / E2) - E1 / exp(E0 / E1) + E0 * gamma(0, E0 / E1) - E0 * gamma(0, E0 / E2)
end


""" 
Electron density after the front using the Lagarkov, Ebert expression with 
Townsend coeff and neglecting the field behind the front.
"""
function ne_integral(Eh)
    α0 = 4.332e5
    E0 = 2e7

    # Dimless correction factor
    w = 1.0
    w * (α0 * co.ϵ0 / co.e) * I0(E0, 0, Eh)
end

    
end # module
