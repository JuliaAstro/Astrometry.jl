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
    
    # Fundamental (Delaunay) arguments
    @inbounds begin
        arg1 = deg2rad(rem(Polynomial(l0_2003A...)(Δt), ARCSECPER2PI) / 3600)
        arg2 = deg2rad(rem(Polynomial(l1_2000A...)(Δt), ARCSECPER2PI) / 3600)
        arg3 = deg2rad(rem(Polynomial(F_2003A...)(Δt), ARCSECPER2PI) / 3600)
        arg4 = deg2rad(rem(Polynomial(D_2000A...)(Δt), ARCSECPER2PI) / 3600)
        arg5 = deg2rad(rem(Polynomial(Ω_2003A...)(Δt), ARCSECPER2PI) / 3600)
    end
    
    ln = ln_2000A_nutation
    la = la_2000A_nutation
    ψl = ϵl = zero(eltype(la))
    @inbounds for i in axes(ln, 1)
        angle = ln[i,1] * arg1 + ln[i,2] * arg2 + ln[i,3] * arg3 + ln[i,4] * arg4 + ln[i,5] * arg5
        angle = rem2pi(angle, RoundToZero)

        s = sin(angle)
        c = cos(angle)
        ψl += (la[i,1] + la[i,2] * Δt) * s + la[i,3] * c
        ϵl += la[i,6] * s + (la[i,4] + la[i,5] * Δt) * c
    end
    
    ###  Planetary Nutation
    
    planet_args = SVector(
        rem2pi(Polynomial(l0_2000A_planet...)(Δt), RoundToZero),
        rem2pi(Polynomial(F_2000A_planet...)(Δt), RoundToZero),
        rem2pi(Polynomial(D_2000A_planet...)(Δt), RoundToZero),
        rem2pi(Polynomial(Ω_2000A_planet...)(Δt), RoundToZero),
        rem2pi(Polynomial(lme_2003...)(Δt), RoundToZero),
        rem2pi(Polynomial(lve_2003...)(Δt), RoundToZero),
        rem2pi(Polynomial(lea_2003...)(Δt), RoundToZero),
        rem2pi(Polynomial(lma_2003...)(Δt), RoundToZero),
        rem2pi(Polynomial(lju_2003...)(Δt), RoundToZero),
        rem2pi(Polynomial(lsa_2003...)(Δt), RoundToZero),
        rem2pi(Polynomial(lur_2003...)(Δt), RoundToZero),
        rem2pi(Polynomial(lne_2003mhb...)(Δt), RoundToZero),
        Polynomial(lge_2003...)(Δt)
    )
    
    pn = pn_2000A_nutation
    pa = pa_2000A_nutation
    ψp = ϵp = zero(eltype(pa))
    @inbounds for i in axes(pn, 1)
        angle = zero(eltype(pn))
        for j in 1:13
            angle += pn[i,j] * planet_args[j]
        end
        angle = rem2pi(angle, RoundToZero)

        s = sin(angle)
        c = cos(angle)
        ψp += pa[i,1] * s + pa[i,2] * c
        ϵp += pa[i,3] * s + pa[i,4] * c
    end
    
    # Convert from 0.1 μas to radians
    FACTOR = deg2rad(1 / 3.6e10)
    ((ψl + ψp) * FACTOR, (ϵl + ϵp) * FACTOR)
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
    #
    # Fundamental (Delaunay) arguments from Simon et al. (1994)
    DEG_ARCSEC = deg2rad(1/3600)
    @inbounds begin
        arg1 = DEG_ARCSEC * rem(Polynomial(l0_2000B...)(Δt), ARCSECPER2PI)
        arg2 = DEG_ARCSEC * rem(Polynomial(l1_2000B...)(Δt), ARCSECPER2PI)
        arg3 = DEG_ARCSEC * rem(Polynomial(F_2000B...)(Δt), ARCSECPER2PI)
        arg4 = DEG_ARCSEC * rem(Polynomial(D_2000B...)(Δt), ARCSECPER2PI)
        arg5 = DEG_ARCSEC * rem(Polynomial(Ω_2000B...)(Δt), ARCSECPER2PI)
    end

    ln = ln_2000B_nutation
    la = la_2000B_nutation
    ψl = ϵl = zero(eltype(la))
    @inbounds for i in axes(ln, 1)
        angle = ln[i,1] * arg1 + ln[i,2] * arg2 + ln[i,3] * arg3 + ln[i,4] * arg4 + ln[i,5] * arg5
        angle = rem2pi(angle, RoundToZero)

        s = sin(angle)
        c = cos(angle)
        ψl += (la[i,1] + la[i,2] * Δt) * s + la[i,3] * c
        ϵl += la[i,6] * s + (la[i,4] + la[i,5] * Δt) * c
    end

    # Convert from 0.1 μas to radians and add planetary correction
    FACTOR = deg2rad(1 / 3.6e10)
    ((ψl + 1e4 * ψ_2000B_planet) * FACTOR, (ϵl + 1e4 * ϵ_2000B_planet) * FACTOR)
end

function iau_2000b_xys(date)
end
