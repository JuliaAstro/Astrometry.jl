####    Astronomy / Ephemerides    ####

"""
    epv00(day1::AbstractFloat, day2::AbstractFloat)

Earth position and velocity, heliocentric and barycentric, with
respect to the Barycentric Celestial Reference System.

# Input

 - `day1`  -- TDB date (Note 1)
 - `day2`  -- TDB date (Note 1)

# Output

 - `pvh`   -- heliocentric Earth position/velocity
 - `pvb`   -- barycentric Earth position/velocity

# Note

1) The TDB date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TDB)=2450123.7 could be expressed in any of these ways, among
   others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 method is best matched to the way the
   argument is handled internally and will deliver the optimum
   resolution.  The MJD method and the date & time methods are both
   good compromises between resolution and convenience.  However, the
   accuracy of the result is more likely to be limited by the
   algorithm itself than the way the date has been expressed.

   n.b. TT can be used instead of TDB in most applications.

2) On return, the arrays pvh and pvb contain the following:

      pvh[0][0]  x       }
      pvh[0][1]  y       } heliocentric position, au
      pvh[0][2]  z       }

      pvh[1][0]  xdot    }
      pvh[1][1]  ydot    } heliocentric velocity, au/d
      pvh[1][2]  zdot    }

      pvb[0][0]  x       }
      pvb[0][1]  y       } barycentric position, au
      pvb[0][2]  z       }

      pvb[1][0]  xdot    }
      pvb[1][1]  ydot    } barycentric velocity, au/d
      pvb[1][2]  zdot    }

   The vectors are with respect to the Barycentric Celestial Reference
   System.  The time unit is one day in TDB.

3) The function is a SIMPLIFIED SOLUTION from the planetary theory
   VSOP2000 (X. Moisson, P. Bretagnon, 2001, Celes. Mechanics &
   Dyn. Astron., 80, 3/4, 205-213) and is an adaptation of original
   Fortran code supplied by P. Bretagnon (private comm., 2000).

4) Comparisons over the time span 1900-2100 with this simplified
   solution and the JPL DE405 ephemeris give the following results:

                              RMS    max
         Heliocentric:
            position error    3.7   11.2   km
            velocity error    1.4    5.0   mm/s

         Barycentric:
            position error    4.6   13.4   km
            velocity error    1.4    4.9   mm/s

   Comparisons with the JPL DE406 ephemeris show that by 1800 and 2200
   the position errors are approximately double their 1900-2100 size.
   By 1500 and 2500 the deterioration is a factor of 10 and by 1000
   and 3000 a factor of 60.  The velocity accuracy falls off at about
   half that rate.

5) It is permissible to use the same array for pvh and pvb, which will
   receive the barycentric values.
"""
function epv00(day1::AbstractFloat, day2::AbstractFloat)
    Δt = ((day1 - JD2000) + day2)/DAYPERYEAR
    @assert abs(Δt) <= 100.0 "Julian day is not between 1990 and 2100."

    # Sun to Earth ecliptic vector
    p_heli = iau_2000_bcrs * SVector(
         ephem_position(sun_earth_x_0, sun_earth_x_1, sun_earth_x_2, Δt),
         ephem_position(sun_earth_y_0, sun_earth_y_1, sun_earth_y_2, Δt),
         ephem_position(sun_earth_z_0, sun_earth_z_1, sun_earth_z_2, Δt))
    v_heli = iau_2000_bcrs * SVector(
         ephem_velocity(sun_earth_x_0, sun_earth_x_1, sun_earth_x_2, Δt),
         ephem_velocity(sun_earth_y_0, sun_earth_y_1, sun_earth_y_2, Δt),
         ephem_velocity(sun_earth_z_0, sun_earth_z_1, sun_earth_z_2, Δt))

    # Barycenter to Earth ecliptic vector
    p_bary = p_heli .+ iau_2000_bcrs * SVector(
         ephem_position(bary_sun_x_0, bary_sun_x_1, bary_sun_x_2, Δt),
         ephem_position(bary_sun_y_0, bary_sun_y_1, bary_sun_y_2, Δt),
         ephem_position(bary_sun_z_0, bary_sun_z_1, bary_sun_z_2, Δt))
    v_bary = v_heli .+ iau_2000_bcrs * SVector(
         ephem_velocity(bary_sun_x_0, bary_sun_x_1, bary_sun_x_2, Δt),
         ephem_velocity(bary_sun_y_0, bary_sun_y_1, bary_sun_y_2, Δt),
         ephem_velocity(bary_sun_z_0, bary_sun_z_1, bary_sun_z_2, Δt))

    (helio = SVector(p_heli, v_heli), bary = SVector(p_bary, v_bary))
end

"""
    moon98(day1::AbstractFloat, day2::AbstractFloat)

Approximate geocentric position and velocity of the Moon.

n.b. Not IAU-endorsed and without canonical status.

# Input

 - `day1`  -- TT date part A (Notes 1,4)
 - `day2`  -- TT date part B (Notes 1,4)

# Output

 - `pv`    -- Moon p,v, GCRS (AU, AU/d, Note 5)

# Note

1) The TT date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TT)=2450123.7 could be expressed in any of these ways, among
   others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 method is best matched to the way the
   argument is handled internally and will deliver the optimum
   resolution.  The MJD method and the date & time methods are both
   good compromises between resolution and convenience.  The limited
   accuracy of the present algorithm is such that any of the methods
   is satisfactory.

2) This function is a full implementation of the algorithm published
   by Meeus (see reference) except that the light-time correction to
   the Moon's mean longitude has been omitted.

3) Comparisons with ELP/MPP02 over the interval 1950-2100 gave RMS
   errors of 2.9 arcsec in geocentric direction, 6.1 km in position
   and 36 mm/s in velocity.  The worst case errors were 18.3 arcsec in
   geocentric direction, 31.7 km in position and 172 mm/s in velocity.

4) The original algorithm is expressed in terms of "dynamical time",
   which can either be TDB or TT without any significant change in
   accuracy.  UT cannot be used without incurring significant errors
   (30 arcsec in the present era) due to the Moon's 0.5 arcsec/sec
   movement.

5) The result is with respect to the GCRS (the same as J2000.0 mean
   equator and equinox to within 23 mas).

6) Velocity is obtained by a complete analytical differentiation of
   the Meeus model.

7) The Meeus algorithm generates position and velocity in mean
   ecliptic coordinates of date, which the present function then
   rotates into GCRS.  Because the ecliptic system is precessing,
   there is a coupling between this spin (about 1.4 degrees per
   century) and the Moon position that produces a small velocity
   contribution.  In the present function this effect is neglected as
   it corresponds to a maximum difference of less than 3 mm/s and
   increases the RMS error by only 0.4%.

# References

Meeus, J., Astronomical Algorithms, 2nd edition, Willmann-Bell, 1998,
p337.

Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M., Francou,
G. & Laskar, J., Astron.Astrophys., 1994, 282, 663
"""
function moon98(day1::AbstractFloat, day2::AbstractFloat)
    Δt = ((day1 - JD2000) + day2)/(100*DAYPERYEAR)

    #  Arguments (radians) and derivatives (radians per Julian century).
    
    #  Moon's mean longitude.
    λm  = deg2rad(rem(Polynomial(λmoon_1994...)(Δt), 360.0))
    dλm = deg2rad(Polynomial(SVector(1., 2., 3., 4.).*λmoon_1994[2:5]...)(Δt))
    
    #  Moon's mean elongation.
    dm  = deg2rad(rem(Polynomial(dmoon_1998...)(Δt), 360.0))
    ddm = deg2rad(Polynomial(SVector(1., 2., 3., 4.).*dmoon_1998[2:5]...)(Δt))

    #  Sun's mean anomaly.
    ls  = deg2rad(rem(Polynomial(lsun_1998...)(Δt), 360.0))
    dls = deg2rad(Polynomial(SVector(1., 2., 3., 4.).*lsun_1998[2:5]...)(Δt))

    #  Moon's mean anomaly.
    lm  = deg2rad(rem(Polynomial(lmoon_1998...)(Δt), 360.0))
    dlm = deg2rad(Polynomial(SVector(1., 2., 3., 4.).*lmoon_1998[2:5]...)(Δt))

    #  Mean distance of the Moon from its ascending node.
    fm  = deg2rad(rem(Polynomial(fmoon_1998...)(Δt), 360.0))
    dfm = deg2rad(Polynomial(SVector(1., 2., 3., 4.).*fmoon_1998[2:5]...)(Δt))

    #  Meeus further arguments.
    a1, da1 = deg2rad.(SVector(Polynomial(a_1...)(Δt), a_1[2]))
    a2, da2 = deg2rad.(SVector(Polynomial(a_2...)(Δt), a_2[2]))
    a3, da3 = deg2rad.(SVector(Polynomial(a_3...)(Δt), a_3[2]))

    #  E-factor
    e  = Polynomial(efac...)(Δt)
    de = Polynomial(SVector(1., 2.).*efac[2:3]...)(Δt)

    #  Arange matrices for vector operations.
    ln  = lr_1998_n
    la  = lr_1998_a
    lre  = [abs(m) == 2 ? e*e    : (abs(m) == 1 ? e  : 1.) for m in ln[:,2]]
    dlre = [abs(m) == 2 ? 2*e*de : (abs(m) == 1 ? de : 0.) for m in ln[:,2]]

    bn   = b_1998_n
    ba   = b_1998_a
    bne  = SVector((abs(m) == 2 ? e*e    : (abs(m) == 1 ? e  : 1.) for m in bn[:,2])...)
    dbne = SVector((abs(m) == 2 ? 2*e*de : (abs(m) == 1 ? de : 0.) for m in bn[:,2])...)

    ϕ, dϕ = ln*SVector(dm, ls, lm, fm), ln*SVector(ddm, dls, dlm, dfm)
    ψ, dψ = SVector(a1, λm-fm, a2), SVector(da1, dλm-dfm, da2)
    ζ, dζ  = bn*SVector(dm, ls, lm, fm), bn*SVector(ddm, dls, dlm, dfm)
    η     = SVector(λm, a3, a1-fm, a1+fm, λm-lm, λm+lm)
    dη    = SVector(dλm, da3, da1-dfm, da1+dfm, dλm-dlm, dλm+dlm)

    #  Longitude, latitude, and distance plus derivatives
    λ  =  λm  + deg2rad(sum(a_l.*sin.(ψ)) + sum(la[:,1].*lre.*sin.(ϕ)))
    dλ = (dλm + deg2rad(sum(a_l.*dψ.*cos.(ψ)) +
                        sum(la[:,1].*(dlre.*sin.(ϕ) .+ lre.*dϕ.*cos.(ϕ))))) /
                        (100*DAYPERYEAR)
    b  = deg2rad(sum(a_b.*sin.(η)) + sum(ba[:,1].*bne.*sin.(ζ)))
    db = deg2rad(sum(a_b.*dη.*cos.(η)) +
                 sum(ba[:,1].*(dbne.*sin.(ζ) .+ bne.*dζ.*cos.(ζ))))/(100*DAYPERYEAR)
    r  = (r0 + sum(la[:,2].*lre.*cos.(ϕ)))/ASTRUNIT
    dr = sum(la[:,2].*(dlre.*cos.(ϕ) .- lre.*dϕ.*sin.(ϕ)))/(ASTRUNIT*100*DAYPERYEAR)

    #  IAU 2006 Fukashima-Williams bias+precession angles
    γB, ϕB, ψB, ϵA = pfw06(day1, day2)
    #  Rotate the Moon position and velocity into GCRS (Note 6).
    R, pv = Rz(-γB)*Rx(-ϕB)*Rz(ψB), s2pv(λ, b, r, dλ, db, dr)
    SVector(R*pv[1], R*pv[2])
end

const KMAX = 10

"""
    plan94(day1::AbstractFloat, day2::AbstractFloat, planet::Integer)

Approximate heliocentric position and velocity of a nominated major
planet: Mercury, Venus, EMB, Mars, Jupiter, Saturn, Uranus or Neptune
(but not the Earth itself).

n.b. Not IAU-endorsed and without canonical status.

# Input

 - `day1`   -- TDB date part A (Note 1)
 - `day2`   -- TDB date part B (Note 1)
 - `planet` -- planet (1=Mercury, 2=Venus, 3=EMB, 4=Mars, 5=Jupiter, 6=Saturn,
                       7=Uranus, 8=Neptune)

# Output

 - `pv`     -- planet p,v (heliocentric, J2000.0, au,au/d)

# Note

1) The date date1+date2 is in the TDB time scale (in practice TT can
   be used) and is a Julian Date, apportioned in any convenient way
   between the two arguments.  For example, JD(TDB)=2450123.7 could be
   expressed in any of these ways, among others:

          date1          date2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 method is best matched to the way the
   argument is handled internally and will deliver the optimum
   resolution.  The MJD method and the date & time methods are both
   good compromises between resolution and convenience.  The limited
   accuracy of the present algorithm is such that any of the methods
   is satisfactory.

2) If an np value outside the range 1-8 is supplied, an error status
   (function value -1) is returned and the pv vector set to zeroes.

3) For np=3 the result is for the Earth-Moon Barycenter.  To obtain
   the heliocentric position and velocity of the Earth, use instead
   the ERFA function eraEpv00.

4) On successful return, the array pv contains the following:

      pv[0][0]   x      }
      pv[0][1]   y      } heliocentric position, au
      pv[0][2]   z      }

      pv[1][0]   xdot   }
      pv[1][1]   ydot   } heliocentric velocity, au/d
      pv[1][2]   zdot   }

   The reference frame is equatorial and is with respect to the mean
   equator and equinox of epoch J2000.0.

5) The algorithm is due to J.L. Simon, P. Bretagnon, J. Chapront,
   M. Chapront-Touze, G. Francou and J. Laskar (Bureau des Longitudes,
   Paris, France).  From comparisons with JPL ephemeris DE102, they
   quote the following maximum errors over the interval 1800-2050:

                   L (arcsec)    B (arcsec)      R (km)

      Mercury          4             1             300
      Venus            5             1             800
      EMB              6             1            1000
      Mars            17             1            7700
      Jupiter         71             5           76000
      Saturn          81            13          267000
      Uranus          86             7          712000
      Neptune         11             1          253000

   Over the interval 1000-3000, they report that the accuracy is no
   worse than 1.5 times that over 1800-2050.  Outside 1000-3000 the
   accuracy declines.

   Comparisons of the present function with the JPL DE200 ephemeris
   give the following RMS errors over the interval 1960-2025:

                    position (km)     velocity (m/s)

      Mercury            334               0.437
      Venus             1060               0.855
      EMB               2010               0.815
      Mars              7690               1.98
      Jupiter          71700               7.70
      Saturn          199000              19.4
      Uranus          564000              16.4
      Neptune         158000              14.4

   Comparisons against DE200 over the interval 1800-2100 gave the
   following maximum absolute differences.  (The results using DE406
   were essentially the same.)

                 L (arcsec)   B (arcsec)     R (km)   Rdot (m/s)

      Mercury        7            1            500       0.7
      Venus          7            1           1100       0.9
      EMB            9            1           1300       1.0
      Mars          26            1           9000       2.5
      Jupiter       78            6          82000       8.2
      Saturn        87           14         263000      24.6
      Uranus        86            7         661000      27.4
      Neptune       11            2         248000      21.4

6) The present ERFA re-implementation of the original Simon et al.
   Fortran code differs from the original in the following respects:

     *  C instead of Fortran.

     *  The date is supplied in two parts.

     *  The result is returned only in equatorial Cartesian form;
        the ecliptic longitude, latitude and radius vector are not
        returned.

     *  The result is in the J2000.0 equatorial frame, not ecliptic.

     *  More is done in-line: there are fewer calls to subroutines.

     *  Different error/warning status values are used.

     *  A different Kepler's-equation-solver is used (avoiding
        use of double precision complex).

     *  Polynomials in t are nested to minimize rounding errors.

     *  Explicit double constants are used to avoid mixed-mode
        expressions.

   None of the above changes affects the result significantly.

7) The returned status indicates the most serious condition
   encountered during execution of the function.  Illegal np is
   considered the most serious, overriding failure to converge, which
   in turn takes precedence over the remote date warning.

# References

Simon, J.L, Bretagnon, P., Chapront, J., Chapront-Touze, M., Francou,
G., and Laskar, J., Astron.Astrophys., 282, 663 (1994).
"""
function plan94(day1::AbstractFloat, day2::AbstractFloat, planet::Integer)
    @assert 1 <= planet <= 9 "Invalid planet: valid range is [1-9]"
    Δt = ((day1 - JD2000) + day2)/(1000*DAYPERYEAR)
    @assert abs(Δt) <= 1.0 "Invalid date: valid range is +/-1000 years"

    #  Compute the mean elements.
    μ = 0.35953620*Δt

    a = Polynomial(a_1994[planet,:]...)(Δt) + 1e-7 * 
        (sum(a_cos_1994[planet,1:8].*cos.(p_1994[planet,1:8].*μ) .+
             a_sin_1994[planet,1:8].*sin.(p_1994[planet,1:8].*μ)) + 
         (a_cos_1994[planet,9]*cos(p_1994[planet,9]*μ) +
          a_sin_1994[planet,9]*sin(p_1994[planet,9]*μ))*Δt)
    λ = deg2rad(Polynomial(((3600., 1., 1.).*λ_1994[planet,:])...)(Δt)/3600.) + 1e-7 *
        (sum(λ_cos_1994[planet,1:8].*cos.(q_1994[planet,1:8].*μ) .+
             λ_sin_1994[planet,1:8].*sin.(q_1994[planet,1:8].*μ)) +
         sum(λ_cos_1994[planet,9:10].*cos.(q_1994[planet,9:10].*μ) .+
             λ_sin_1994[planet,9:10].*sin.(q_1994[planet,9:10].*μ))*Δt)
    e = Polynomial(e_1994[planet,:]...)(Δt)
    p = rem2pi(deg2rad(Polynomial(((3600., 1., 1.).*π_1994[planet,:])...)(Δt)/3600.),
               RoundToZero)
    i = deg2rad(Polynomial(((3600., 1., 1.).*i_1994[planet,:])...)(Δt)/3600.)
    ω = rem2pi(deg2rad(Polynomial(((3600., 1., 1.).*ω_1994[planet,:])...)(Δt)/3600.),
               RoundToZero)

    #  Iterative solution of Kepler's equation to get eccentric anomaly.
    am = rem2pi(λ, RoundToZero) - p
    ae, dae = am + e*sin(am), 1.0
    for k=1:KMAX
        dae = (am - ae + e*sin(ae))/(1 - e*cos(ae))
        ae += dae
        if abs(dae) <= 1e-12
            break
        elseif k == KMAX
            @warn "Maximum iterations reached."
        end
    end

    #  True anomaly
    at = 2.0*atan(sqrt((1.0+e)/(1.0-e))*sin(ae/2.), cos(ae/2.))

    #  Distance (AU) and speed (radians/day).
    r, v = a*(1 - e*cos(ae)), GK*sqrt((1 + 1/mass_1994[planet])/a^3)
    xm2 = 2.0*sin(i/2.0)*(sin(ω)*cos(at + p) - cos(ω)*sin(at + p))
    xms, xmc = a/sqrt(1-e*e).*(e*sin(p) + sin(at + p), e*cos(p) + cos(at + p))

    #  Position (J2000.0 ecliptic x, y, z in AU).
    px = r*(cos(at + p) - xm2*sin(i/2.0)*sin(ω))
    py = r*(sin(at + p) + xm2*sin(i/2.0)*cos(ω))
    pz = -r*xm2*cos(i/2.0)

    #  Velocity (J2000.0 ecliptic v_x, v_y, v_z in AU/day).
    vx = -v*((1 - 2*(sin(i/2)*sin(ω))^2)*xms - 2*sin(i/2)^2*sin(ω)*cos(ω)*xmc)
    vy =  v*((1 - 2*(sin(i/2)*cos(ω))^2)*xmc - 2*sin(i/2)^2*sin(ω)*cos(ω)*xms)
    vz =  v*2*cos(i/2)*sin(i/2)*(sin(ω)*xms + cos(ω)*xmc)

    #  Rotate to equatorial.
    SVector(SVector(px, COSEPS*py - SINEPS*pz, SINEPS*py + COSEPS*pz),
            SVector(vx, COSEPS*vy - SINEPS*vz, SINEPS*vy + COSEPS*vz))
end
