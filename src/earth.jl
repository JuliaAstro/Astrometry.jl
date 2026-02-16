"""
    equinox(date, model=:IAU2006A)
"""
function equinox(date, model=:IAU2006A)
    models = [:IAU2000A, :IAU2000B, :IAU2006]
    @assert model in models "Incorrect model type"
    
    return if model == :IAU2000A
        ϵA = iau_1980_obliquity(date) + ϵ_2000_corr*(date - JD2000)/(100*DAYPERYEAR)
        ψn = iau_2000a_nutation(date)[1]
        ψn * cos(ϵA) + equinox_complement(date)
    elseif model == :IAU2000B
        ϵA = iau_1980_obliquity(date) + ϵ_2000_corr*(date - JD2000)/(100*DAYPERYEAR)
        ψn = iau_2000b_nutation(date)[1]
        ψn * cos(ϵA) + equinox_complement(date)
    elseif model == :IAU2006
        rem2pi(iau_2006a_gst(0.0, date) - iau_2006_gmst(0.0, date), RoundNearest)
    end
end

"""
    precession_nutation(date; model=:IAU2006A)


"""
function precession_nutation(date; model=:IAU2006A)
    models = [:IAU1980, :IAU2000A, :IAU2000B, :IAU2006A]
    @assert model in models "Incorrect model type"
    
    return if model == :IAU1980
        ζ, z, θ = iau_1976_precession(JD2000, date)
        ϵA = iau_1980_obliquity(date)
        ψ, ϵ = iau_1980_nutation(date)
        Rx(-(ϵA + ϵ)) * Rz(-ψ) * Rx(ϵA) * Rz(-z) * Ry(θ) * Rz(-ζ)
    elseif model == :IAU2000A
        # Frame bias angles: GCRS to 2000.0
        rab, ψb, ϵb, ϵ0 = deg2rad(1/3600).*(icrs_ra_2000, ψ_bias_2000, ϵ_bias_2000, ϵ0_2000)
        ψ, ω, χ = iau_2000_precession(date)
        ϵA = iau_1980_obliquity(date) + deg2rad(1/3600)*ϵ_corr_2000*(date - JD2000)/(100*DAYPERYEAR)
        ψn, ϵn = iau_2000a_nutation(date)
        Rx(-(ϵA + ϵn)) * Rz(-ψn) * Rx(ϵA) * Rz(χ) * Rx(-ω) * Rz(-ψ) * Rx(ϵ0) * Rx(-ϵb) * Ry(ψb*sin(ϵ0)) * Rz(rab)
    elseif model == :IAU2000B
        rab, ψb, ϵb, ϵ0 = deg2rad(1/3600).*(icrs_ra_2000, ψ_bias_2000, ϵ_bias_2000, ϵ0_2000)
        ψ, ω, χ = iau_2000_precession(date)
        ϵA = iau_1980_obliquity(date) + deg2rad(1/3600)*ϵ_corr_2000*(date - JD2000)/(100*DAYPERYEAR)
        ψn, ϵn = iau_2000b_nutation(date)
        rmat = Rx(-(ϵA + ϵn)) * Rz(-ψn) * Rx(ϵA) * Rz(χ) * Rx(-ω) * Rz(-ψ) * Rx(ϵ0) * Rx(-ϵb) * Ry(ψb*sin(ϵ0)) * Rz(rab)
    elseif model == :IAU2006A
        γB, ϕB, ψB, ϵB = iau_2006_precession(date)
        ψn, ϵn = iau_2006a_nutation(date)
        Rx(-(ϵB + ϵn)) * Rz(-(ψB + ψn)) * Rx(ϕB) * Rz(γB)
    end
end
