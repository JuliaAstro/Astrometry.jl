# IAU 2000 Model

include("constants2000.jl")

"""
    iau_2000_earth_rotation(date)

Earth rotation angle (IAU 2000)
"""
function iau_2000_earth_rotation(date)
    
    rem2pi(mod(date) + Polynomial(era_2000...)(date - JD2000))
end

function iau_2000_equinox_complement(date)

    Δt = (date - JD2000)/(100*DAYPERYEAR)

    ϕ = deg2rad.(rem.([
        Polynomial(l0_2003A...)(Δt), Polynomial(l1_2003A...)(Δt),
        Polynomial( F_2003A...)(Δt), Polynomial( D_2000A...)(Δt),
        Polynomial( Ω_2003A...)(Δt)], ARCSECPER2PI)./3600)
    append!(ϕ, [
        Polynomial( lve_2003...)(Δt), Polynomial( lea_2003...)(Δt),
        Polynomial( lge_2003...)(Δt)])
    
    en0 = ϕ0_2000_equinox
    ea0 = a0_2000_equinox
    en1 = ϕ1_2000_equinox
    ea1 = a1_2000_equinox
    deg2rad((sum(ea0[:,1].*sin.(en0*ϕ) .+ ea0[:,2].*cos.(en0*ϕ)) +
             sum(ea1[:,1].*sin.(en1*ϕ) .+ ea1[:,2].*cos.(en1*ϕ))*Δt)/3600)
end

function iau_2000_gmst(ut, tt)
end

function ephem_position(coef0, coef1, coef2, Δt)

    A0, ϕ0, ν0 = [coef0[j,:] for j=1:3]
    A1, ϕ1, ν1 = [coef1[j,:] for j=1:3]
    A2, ϕ2, ν2 = [coef2[j,:] for j=1:3]

    (sum(A0 .* cos.(ϕ0 .+ ν0 .* Δt)) +
     sum(A1 .* cos.(ϕ1 .+ ν1 .* Δt))*Δt +
     sum(A2 .* cos.(ϕ2 .+ ν2 .* Δt))*Δt^2)
end

function ephem_velocity(coef0, coef1, coef2, Δt)

    A0, ϕ0, ν0 = coef0[1,:], coef0[2,:], coef0[3,:]
    A1, ϕ1, ν1 = coef1[1,:], coef1[2,:], coef1[3,:]
    A2, ϕ2, ν2 = coef2[1,:], coef2[2,:], coef2[3,:]
    
    (-sum(A0 .*  ν0 .* sin.(ϕ0 .+ ν0 .* Δt)) +
      sum(A1 .* (cos.(ϕ1 .+ ν1 .* Δt) .- ν1 .* Δt .* sin.(ϕ1 .+ ν1 .* Δt))) +
      sum(A2 .* (2 .* cos.(ϕ2 .+ ν2 .* Δt) .-
                 ν2 .* Δt .* sin.(ϕ2 .+ ν2 .* Δt)))*Δt)/DAYPERYEAR
end

include("ephemerisDE405.jl")

"""
    iau_2000_earth_position(date; frame=:barycenter)
"""
function iau_2000_position(date; frame=:barycenter)

    @assert frame in [:barycenter, :heliocenter]
    Δt = (date - JD2000)/DAYPERYEAR
    @assert abs(Δt) <= 100.0 "Julian day is not between 1990 and 2100."

    # Sun to Earth ecliptic vector
    position = [
        ephem_position(sun_earth_x_0, sun_earth_x_1, sun_earth_x_2, Δt),
        ephem_position(sun_earth_y_0, sun_earth_y_1, sun_earth_y_2, Δt),
        ephem_position(sun_earth_z_0, sun_earth_z_1, sun_earth_z_2, Δt)]
    velocity = [
        ephem_velocity(sun_earth_x_0, sun_earth_x_1, sun_earth_x_2, Δt),
        ephem_velocity(sun_earth_y_0, sun_earth_y_1, sun_earth_y_2, Δt),
        ephem_velocity(sun_earth_z_0, sun_earth_z_1, sun_earth_z_2, Δt)]

    # Barycenter to Earth ecliptic vector
    if frame == :barycenter
        position .+= [
            ephem_position(bary_sun_x_0, bary_sun_x_1, bary_sun_x_2, Δt),
            ephem_position(bary_sun_y_0, bary_sun_y_1, bary_sun_y_2, Δt),
            ephem_position(bary_sun_z_0, bary_sun_z_1, bary_sun_z_2, Δt)]
        velocity .+= [
            ephem_velocity(bary_sun_x_0, bary_sun_x_1, bary_sun_x_2, Δt),
            ephem_velocity(bary_sun_y_0, bary_sun_y_1, bary_sun_y_2, Δt),
            ephem_velocity(bary_sun_z_0, bary_sun_z_1, bary_sun_z_2, Δt)]
    end
    (iau_2000_bcrs*position, iau_2000_bcrs*velocity)
end

"""

Frame bias and precession, IAU 2000
"""
function iau_2000_precession(date)

    Δt = (date - JD2000)/(100*DAYPERYEAR)
    
    ψ = Polynomial(ψ_1977...)(Δt) + ψ_corr_2000*Δt
    ϵ = Polynomial(ω_1977...)(Δt) + ϵ_corr_2000*Δt
    χ  = Polynomial(χ_1977...)(Δt)

    deg2rad(1/3600).*(ψ, ϵ, χ)
end

function iau_2000_tio_locator(date)
end

#    IAU 2000A Model

include("constants2000A.jl")
include("constants2000B.jl")

function iau_2000a_cio_locator(date)
end

function iau_2000a_crs_cis(date)
end

function iau_2000a_crs_trs(date)
end

function iau_2000a_gst(ut, tt)
end

function iau_2000a_nutation(date)

    Δt = (date - JD2000)/(100*DAYPERYEAR)

    ###  Luni-solar Nutation
    
    # Fundamental (Delauney) arguments
    ϕ = deg2rad.(rem.(
        [Polynomial(l0_2003A...)(Δt), Polynomial(l1_2000A...)(Δt),
         Polynomial( F_2003A...)(Δt), Polynomial( D_2000A...)(Δt),
         Polynomial( Ω_2003A...)(Δt)], ARCSECPER2PI)/3600)

    lϕ  = rem2pi.(vcat([t.n' for t in iau_2000A_nutation_lunisolar_series]...)*ϕ, RoundToZero)
    la  = vcat([t.a' for t in iau_2000A_nutation_lunisolar_series]...)
    ψl = sum((la[:,1] .+ la[:,2].*Δt).*sin.(lϕ) .+ la[:,3].*cos.(lϕ))
    ϵl = sum(la[:,6].*sin.(lϕ) .+ (la[:,4] .+ la[:,5].*Δt).*cos.(lϕ))

    ###  Planetary Nutation

    ϕ = rem2pi.([
        Polynomial(l0_2000A_planet...)(Δt), Polynomial( F_2000A_planet...)(Δt),
        Polynomial( D_2000A_planet...)(Δt), Polynomial( Ω_2000A_planet...)(Δt),
        Polynomial(lme_2003...)(Δt), Polynomial(lve_2003...)(Δt),
        Polynomial(lea_2003...)(Δt), Polynomial(lma_2003...)(Δt),
        Polynomial(lju_2003...)(Δt), Polynomial(lsa_2003...)(Δt),
        Polynomial(lur_2003...)(Δt), Polynomial(lne_2003mhb...)(Δt)], RoundToZero)
    push!(ϕ, Polynomial(lge_2003...)(Δt))

    pϕ = rem2pi.(vcat([t.n' for t in iau_2000A_nutation_planetary_series]...)*ϕ, RoundToZero)
    pa = vcat([t.a' for t in iau_2000A_nutation_planetary_series]...)
    ψp = sum(pa[:,1].*sin.(pϕ) .+ pa[:,2].*cos.(pϕ))
    ϵp = sum(pa[:,3].*sin.(pϕ) .+ pa[:,4].*cos.(pϕ))

    deg2rad.((ψl + ψp, ϵl + ϵp)./3.6e10)
end

function iau_2000a_xys(date)
end

#   IAU 2000B Model

function iau_2000b_cio_locator(date)
end

function iau_2000b_crs_cis(date)
end

function iau_2000b_crs_trs(date)
end

function iau_2000b_gst(ut, tt)
end

function iau_2000b_nutation(date)

    Δt = (date - JD2000)/(100*DAYPERYEAR)

    ###  Luni-solar Nutation
    
    # Fundamental (Delauney) arguments from Simon et al. (1994)
    ϕ = deg2rad(1/3600).*rem.([Polynomial(l0_2000B...)(Δt),
        Polynomial(l1_2000B...)(Δt), Polynomial( F_2000B...)(Δt),
        Polynomial( D_2000B...)(Δt), Polynomial( Ω_2000B...)(Δt)], ARCSECPER2PI)

    lϕ = rem2pi.(vcat([t.n' for t in iau_2000B_nutation_lunisolar_series]...)*ϕ, RoundToZero)
    la = vcat([t.a' for t in iau_2000B_nutation_lunisolar_series]...)
    ψl = sum((la[:,1] .+ la[:,2].*Δt).*sin.(lϕ) .+ la[:,3].*cos.(lϕ))
    ϵl = sum(la[:,6].*sin.(lϕ) .+ (la[:,4] .+ la[:,5].*Δt).*cos.(lϕ))

    deg2rad.((ψl + 1e4*ψ_2000B_planet, ϵl + 1e4*ϵ_2000B_planet)./3.6e10)
end

function iau_2000b_xys(date)
end
