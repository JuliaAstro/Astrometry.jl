#### Astronomy / Rotation and Time

"""
    ee00(day1::AbstractFloat, day2::AbstractFloat, ϵA::AbstractFloat, ψ::AbstractFloat)

The equation of the equinoxes, compatible with IAU 2000 resolutions,
given the nutation in longitude and the mean obliquity.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date
 - `epsa`  -- mean obliquity (Note 2)
 - `dpsi`  -- nutation in longitude (Note 3)

# Output

 - `ee`    --  equation of the equinoxes (Note 4)

# Note

1) The TT date day1+day2 is a Julian Date, apportioned in any
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
   good compromises between resolution and convenience.

2) The obliquity, in radians, is mean of date.

3) The result, which is in radians, operates in the following sense:

      Greenwich apparent ST = GMST + equation of the equinoxes

4) The result is compatible with the IAU 2000 resolutions.  For
   further details, see IERS Conventions 2003 and Capitaine et al.
   (2002).

# References

Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
implement the IAU 2000 definition of UT1", Astronomy & Astrophysics,
406, 1135-1149 (2003)

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function ee00(day1::AbstractFloat, day2::AbstractFloat, ϵA::AbstractFloat, ψ::AbstractFloat)
    #  Equation of the equinoxes
    ψ*cos(ϵA) + eect00(day1, day2)
end

"""
    ee00a(day1::AbstractFloat, day2::AbstractFloat)

Equation of the equinoxes, compatible with IAU 2000 resolutions.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `ee`    -- equation of the equinoxes (Note 2)

# Note

1) The TT date day1+day2 is a Julian Date, apportioned in any
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
   good compromises between resolution and convenience.

2) The result, which is in radians, operates in the following sense:

      Greenwich apparent ST = GMST + equation of the equinoxes

3) The result is compatible with the IAU 2000 resolutions.  For
   further details, see IERS Conventions 2003 and Capitaine et al.
   (2002).

# References

Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
implement the IAU 2000 definition of UT1", Astronomy & Astrophysics,
406, 1135-1149 (2003).

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004).
"""
function ee00a(day1::AbstractFloat, day2::AbstractFloat)
    #  IAU 2000 precession rate adjustment, mean obliquity, consistent with
    #  IAU 2000 precession-nutation, nutation in longitude, and solve for
    #  equation of the equinoxes.
    ee00(day1, day2, obl80(day1, day2)+pr00(day1, day2)[:ϵ],
         nut00a(day1, day2)[:ψ])
end

"""
    ee00b(day1::AbstractFloat, day2::AbstractFloat)

Equation of the equinoxes, compatible with IAU 2000 resolutions but
using the truncated nutation model IAU 2000B.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `ee`    -- equation of the equinoxes (Note 2)

# Note

1) The TT date day1+day2 is a Julian Date, apportioned in any
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
   good compromises between resolution and convenience.

2) The result, which is in radians, operates in the following sense:

      Greenwich apparent ST = GMST + equation of the equinoxes

3) The result is compatible with the IAU 2000 resolutions except that
   accuracy has been compromised (1 mas) for the sake of speed.  For
   further details, see McCarthy & Luzum (2003), IERS Conventions 2003
   and Capitaine et al. (2003).

# References

Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
implement the IAU 2000 definition of UT1", Astronomy & Astrophysics,
406, 1135-1149 (2003)

McCarthy, D.D. & Luzum, B.J., "An abridged model of the
precession-nutation of the celestial pole", Celestial Mechanics &
Dynamical Astronomy, 85, 37-49 (2003)

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function ee00b(day1::AbstractFloat, day2::AbstractFloat)
    #  IAU 2000 precession rate adjustment, mean obliquity, consistent with
    #  IAU 2000 precession-nutation, nutation in longitude, and solve for
    #  equation of the equinoxes.
    ee00(day1, day2, obl80(day1, day2)+pr00(day1, day2)[:ϵ],
         nut00b(day1, day2)[:ψ])
end

"""
    ee06a(day1::AbstractFloat, day2::AbstractFloat)

Equation of the equinoxes, compatible with IAU 2000 resolutions and
IAU 2006/2000A precession-nutation.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `ee`    -- equation of the equinoxes (Note 2)

# Note

1) The TT date day1+day2 is a Julian Date, apportioned in any
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
   good compromises between resolution and convenience.

2) The result, which is in radians, operates in the following sense:

      Greenwich apparent ST = GMST + equation of the equinoxes

# References

McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003), IERS
Technical Note No. 32, BKG
"""
function ee06a(day1::F, day2::F) where F<:AbstractFloat
    rem2pi(gst06a(0., 0., day1, day2) - gmst06(0., 0., day1, day2), RoundNearest)
end

"""
    eect00(day1::AbstractFloat, day2::AbstractFloat)

Equation of the equinoxes complementary terms, consistent with IAU
2000 resolutions.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `eect`  -- complementary terms (Note 2)

# Note

1) The TT date day1+day2 is a Julian Date, apportioned in any
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
   good compromises between resolution and convenience.

2) The "complementary terms" are part of the equation of the equinoxes
   (EE), classically the difference between apparent and mean Sidereal
   Time:

      GAST = GMST + EE

   with:

      EE = dpsi * cos(eps)

   where dpsi is the nutation in longitude and eps is the obliquity of
   date.  However, if the rotation of the Earth were constant in an
   inertial frame the classical formulation would lead to apparent
   irregularities in the UT1 timescale traceable to side- effects of
   precession-nutation.  In order to eliminate these effects from UT1,
   "complementary terms" were introduced in 1994 (IAU, 1994) and took
   effect from 1997 (Capitaine and Gontier, 1993):

      GAST = GMST + CT + EE

   By convention, the complementary terms are included as part of the
   equation of the equinoxes rather than as part of the mean Sidereal
   Time.  This slightly compromises the "geometrical" interpretation
   of mean sidereal time but is otherwise inconsequential.

   The present function computes CT in the above expression,
   compatible with IAU 2000 resolutions (Capitaine et al., 2002, and
   IERS Conventions 2003).

# References

Capitaine, N. & Gontier, A.-M., Astron.Astrophys., 275, 645-650 (1993)

Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
implement the IAU 2000 definition of UT1", Astron.Astrophys., 406,
1135-1149 (2003)

IAU Resolution C7, Recommendation 3 (1994)

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function eect00(day1::F, day2::F) where F<:AbstractFloat
    #  Interval between fundamental epoch J2000.0 and current date.
    Δt = ((day1 - JD2000) + day2)/(100*DAYPERYEAR)

    #  Fundamental arguments (from IERS Conventions 2003).
    ϕ = SVector(fal03(Δt), falp03(Δt), faf03(Δt), fad03(Δt), faom03(Δt), fave03(Δt), fae03(Δt), fapa03(Δt))

    #  Evaluate the EE complementary terms.
    a0 = a0_2000_equinox
    a1 = a1_2000_equinox

    sum0 = sum1 = zero(eltype(a0))
    @inbounds for i in axes(a0, 1)
        angle = zero(eltype(ϕ0_2000_equinox))
        for j in axes(ϕ0_2000_equinox, 2)
            angle += ϕ0_2000_equinox[i,j] * ϕ[j]
        end
        sum0 += a0[i,1] * sin(angle) + a0[i,2] * cos(angle)
    end

    @inbounds for i in axes(a1, 1)
        angle = zero(eltype(ϕ1_2000_equinox))
        for j in axes(ϕ1_2000_equinox, 2)
            angle += ϕ1_2000_equinox[i,j] * ϕ[j]
        end
        sum1 += a1[i,1] * sin(angle) + a1[i,2] * cos(angle)
    end

    deg2rad((sum0 + sum1 * Δt) / 3600.0)
end

"""
    eqeq94(day1::AbstractFloat, day2::AbstractFloat)

Equation of the equinoxes, IAU 1994 model.

# Input

 - `day1`  -- TDB date (Note 1)
 - `day2`  -- TDB date

# Output

 - `ee`    -- equation of the equinoxes (Note 2)

# Note

1) The date day1+day2 is a Julian Date, apportioned in any convenient
   way between the two arguments.  For example, JD(TT)=2450123.7 could
   be expressed in any of these ways, among others:

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
   good compromises between resolution and convenience.

2) The result, which is in radians, operates in the following sense:

      Greenwich apparent ST = GMST + equation of the equinoxes

# References

IAU Resolution C7, Recommendation 3 (1994).

Capitaine, N. & Gontier, A.-M., 1993, Astron.Astrophys., 275, 645-650.
"""
function eqeq94(day1::AbstractFloat, day2::AbstractFloat)
    #  Interval between fundamental epoch J2000.0 and given date
    #  (Julian centuries).
    Δt = ((day1 - JD2000) + day2)/(100*DAYPERYEAR)
    #  Longitude of the mean ascending node of the lunar orbit on the
    #  ecliptic, measured from the mean equinox of date.
    ω = rem2pi(deg2rad(Polynomial(l_1994...)(Δt)/3600) + 2π*rem(-5.0*Δt, 1.),
               RoundNearest)
    #  Nutation components, mean obliquity, and equation of the equinoxes
    nut80(day1, day2)[:ψ]*cos(obl80(day1, day2)) +
        deg2rad(sum(equinox_1994.*(sin(ω), sin(2ω)))/3600.0)
end

"""
    era00(day1::AbstractFloat, day2::AbstractFloat)

Earth rotation angle (IAU 2000 model).

# Input

 - `day1`  -- UT1 as Julian Date (see note)
 - `day2`  -- ... as Julian Date

# Output

 - `era`   -- Earth rotation angle (radians), range 0-2pi

# Note

1) The UT1 date dj1+dj2 is a Julian Date, apportioned in any
   convenient way between the arguments dj1 and dj2.  For example,
   JD(UT1)=2450123.7 could be expressed in any of these ways, among
   others:

           dj1            dj2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 and MJD methods are good compromises between
   resolution and convenience.  The date & time method is best matched
   to the algorithm used: maximum precision is delivered when the dj1
   argument is for 0hrs UT1 on the day in question and the dj2
   argument lies in the range 0 to 1, or vice versa.

2) The algorithm is adapted from Expression 22 of Capitaine et al.
   2000.  The time argument has been expressed in days directly, and,
   to retain precision, integer contributions have been eliminated.
   The same formulation is given in IERS Conventions (2003), Chap. 5,
   Eq. 14.

# References

Capitaine N., Guinot B. and McCarthy D.D, 2000, Astron.  Astrophys.,
355, 398-405.

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function era00(day1::AbstractFloat, day2::AbstractFloat)
    if day1 < day2
        Δt = day1 + (day2 - JD2000)
        fr = rem(day1, 1.0) + rem(day2, 1.0)
    else
        Δt = day2 + (day1 - JD2000)
        fr = rem(day2, 1.0) + mod(day1, 1.0)
    end
    mod2pi(2π*(fr + Polynomial(era_2000...)(Δt)))
end

"""
    gmst00(ut1::AbstractFloat, ut2::AbstractFloat, tt1::AbstractFloat, tt2::AbstractFloat)

Greenwich mean sidereal time (model consistent with IAU 2000
resolutions).

# Input

 - `uta`   -- UT1 as Julian Date (Notes 1,2)
 - `utb`   -- ... as Julian Date
 - `tta`   -- TT as Julian Date (Notes 1,2)
 - `ttb`   -- ... as Julian Date

# Output

 - `gmst`  -- Greenwich mean sidereal time (radians)

# Note

1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
   Julian Dates, apportioned in any convenient way between the
   argument pairs.  For example, JD(UT1)=2450123.7 could be expressed
   in any of these ways, among others:

          Part A         Part B

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable (in the case of UT; the TT is not at all critical in
   this respect).  The J2000 and MJD methods are good compromises
   between resolution and convenience.  For UT, the date & time method
   is best matched to the algorithm that is used by the Earth Rotation
   Angle function, called internally: maximum precision is delivered
   when the uta argument is for 0hrs UT1 on the day in question and
   the utb argument lies in the range 0 to 1, or vice versa.

2) Both UT1 and TT are required, UT1 to predict the Earth rotation and
   TT to predict the effects of precession.  If UT1 is used for both
   purposes, errors of order 100 microarcseconds result.

3) This GMST is compatible with the IAU 2000 resolutions and must be
   used only in conjunction with other IAU 2000 compatible components
   such as precession-nutation and equation of the equinoxes.

4) The result is returned in the range 0 to 2pi.

5) The algorithm is from Capitaine et al. (2003) and IERS Conventions
   2003.

# References

Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
implement the IAU 2000 definition of UT1", Astronomy & Astrophysics,
406, 1135-1149 (2003)

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function gmst00(ut1::AbstractFloat, ut2::AbstractFloat, tt1::AbstractFloat, tt2::AbstractFloat)
    #  TT Julian centuries since J2000.0.
    mod2pi(era00(ut1, ut2) + deg2rad(Polynomial(gmst_2000...)(
        ((tt1 - JD2000) + tt2)/(100*DAYPERYEAR))/3600.0))
end

"""
    gmst06(ut1::AbstractFloat, ut2::AbstractFloat, tt1::AbstractFloat, tt2::AbstractFloat)

Greenwich mean sidereal time (consistent with IAU 2006 precession).

# Input

 - `uta`   -- UT1 as Julian Date (Notes 1,2)
 - `utb`   -- ... as Julian Date
 - `tta`   -- TT as Julian Date (Notes 1,2)
 - `ttb`   -- ... as Julian Date

# Output

 - `gmst`  -- Greenwich mean sidereal time (radians)

# Note

1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
   Julian Dates, apportioned in any convenient way between the
   argument pairs.  For example, JD=2450123.7 could be expressed in
   any of these ways, among others:

          Part A        Part B

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable (in the case of UT; the TT is not at all critical in
   this respect).  The J2000 and MJD methods are good compromises
   between resolution and convenience.  For UT, the date & time method
   is best matched to the algorithm that is used by the Earth rotation
   angle function, called internally: maximum precision is delivered
   when the uta argument is for 0hrs UT1 on the day in question and
   the utb argument lies in the range 0 to 1, or vice versa.

2) Both UT1 and TT are required, UT1 to predict the Earth rotation and
   TT to predict the effects of precession.  If UT1 is used for both
   purposes, errors of order 100 microarcseconds result.

3) This GMST is compatible with the IAU 2006 precession and must not
   be used with other precession models.

4) The result is returned in the range 0 to 2pi.

# References

Capitaine, N., Wallace, P.T. & Chapront, J., 2005,
Astron.Astrophys. 432, 355
"""
function gmst06(ut1::AbstractFloat, ut2::AbstractFloat, tt1::AbstractFloat, tt2::AbstractFloat)
    #  TT Julian centuries since J2000.0
    Δt = ((tt1 - JD2000) + tt2)/(100*DAYPERYEAR)
    #  Greenwich mean sidereal time, IAU 2006
    rem2pi(era00(ut1, ut2) + deg2rad(Polynomial(gmst_2006...)(Δt)/3600),
           RoundNearest)
end

"""
    gmst82(day1::AbstractFloat, day2::AbstractFloat)

Universal Time to Greenwich mean sidereal time (IAU 1982 model).

# Input

 - `day1`   -- UT1 Julian Date (see note)
 - `day2`   -- ... Julian Date

# Output

 - `gmst`   -- Greenwich mean sidereal time (radians)

# Note

1) The UT1 date dj1+dj2 is a Julian Date, apportioned in any
   convenient way between the arguments dj1 and dj2.  For example,
   JD(UT1)=2450123.7 could be expressed in any of these ways, among
   others:

           dj1            dj2

       2450123.7          0          (JD method)
        2451545        -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5         0.2         (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 and MJD methods are good compromises between
   resolution and convenience.  The date & time method is best matched
   to the algorithm used: maximum accuracy (or, at least, minimum
   noise) is delivered when the dj1 argument is for 0hrs UT1 on the
   day in question and the dj2 argument lies in the range 0 to 1, or
   vice versa.

2) The algorithm is based on the IAU 1982 expression.  This is always
   described as giving the GMST at 0 hours UT1.  In fact, it gives the
   difference between the GMST and the UT, the steady
   4-minutes-per-day drawing-ahead of ST with respect to UT.  When
   whole days are ignored, the expression happens to equal the GMST at
   0 hours UT1 each day.

3) In this function, the entire UT1 (the sum of the two arguments dj1
   and dj2) is used directly as the argument for the standard formula,
   the constant term of which is adjusted by 12 hours to take account
   of the noon phasing of Julian Date.  The UT1 is then added, but
   omitting whole days to conserve accuracy.

# References

Transactions of the International Astronomical Union, XVIII B, 67
(1983).

Aoki et al., Astron.Astrophys., 105, 359-361 (1982).
"""
function gmst82(day1::AbstractFloat, day2::AbstractFloat)
    #  Julian centuries since fundamental epoch.
    Δt = (day1 < day2 ? (day2-JD2000)+day1 : (day1-JD2000)+day2)/(100*DAYPERYEAR)
    #  Fractional part of JD(UT1) (seconds) and GMST at this UT1.
    mod2pi(deg2rad((Polynomial(gmst_1982...)(Δt) +
                    SECPERDAY*(rem(day1, 1.0) + rem(day2, 1.0)))/240.0))
end

"""
    gst00a(ut1::AbstractFloat, ut2::AbstractFloat, tt1::AbstractFloat, tt2::AbstractFloat)

Greenwich apparent sidereal time (consistent with IAU 2000
resolutions).

# Input

 - `uta`   -- UT1 as Julian Date (Notes 1,2)
 - `utb`   -- ... as Julian Date
 - `tta`   -- TT as Julian Date (Notes 1,2)
 - `ttb`   -- ... as Julian Date

# Output

 - `gst`   -- Greenwich apparent sidereal time (radians)

# Note

1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
   Julian Dates, apportioned in any convenient way between the
   argument pairs.  For example, JD(UT1)=2450123.7 could be expressed
   in any of these ways, among others:

           uta            utb

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable (in the case of UT; the TT is not at all critical in
   this respect).  The J2000 and MJD methods are good compromises
   between resolution and convenience.  For UT, the date & time method
   is best matched to the algorithm that is used by the Earth Rotation
   Angle function, called internally: maximum precision is delivered
   when the uta argument is for 0hrs UT1 on the day in question and
   the utb argument lies in the range 0 to 1, or vice versa.

2) Both UT1 and TT are required, UT1 to predict the Earth rotation and
   TT to predict the effects of precession-nutation.  If UT1 is used
   for both purposes, errors of order 100 microarcseconds result.

3) This GAST is compatible with the IAU 2000 resolutions and must be
   used only in conjunction with other IAU 2000 compatible components
   such as precession-nutation.

4) The result is returned in the range 0 to 2pi.

5) The algorithm is from Capitaine et al. (2003) and IERS Conventions
   2003.

# References

Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
implement the IAU 2000 definition of UT1", Astronomy & Astrophysics,
406, 1135-1149 (2003)

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function gst00a(ut1::AbstractFloat, ut2::AbstractFloat, tt1::AbstractFloat, tt2::AbstractFloat)
    mod2pi(gmst00(ut1, ut2, tt1, tt2) + ee00a(tt1, tt2))
end

"""
    gst00b(ut1::AbstractFloat, ut2::AbstractFloat)

Greenwich apparent sidereal time (consistent with IAU 2000 resolutions
but using the truncated nutation model IAU 2000B).

# Input

 - `uta`   -- UT1 as Julian Date (Notes 1,2)
 - `utb`   -- ... as Julian Date

# Output

 - `gst`   -- Greenwich apparent sidereal time (radians)

# Note

1) The UT1 date uta+utb is a Julian Date, apportioned in any
   convenient way between the argument pair.  For example,
   JD(UT1)=2450123.7 could be expressed in any of these ways, among
   others:

           uta            utb

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 and MJD methods are good compromises between
   resolution and convenience.  For UT, the date & time method is best
   matched to the algorithm that is used by the Earth Rotation Angle
   function, called internally: maximum precision is delivered when
   the uta argument is for 0hrs UT1 on the day in question and the utb
   argument lies in the range 0 to 1, or vice versa.

2) The result is compatible with the IAU 2000 resolutions, except that
   accuracy has been compromised for the sake of speed and convenience
   in two respects:

   . UT is used instead of TDB (or TT) to compute the precession
     component of GMST and the equation of the equinoxes.  This
     results in errors of order 0.1 mas at present.

   . The IAU 2000B abridged nutation model (McCarthy & Luzum, 2003)
     is used, introducing errors of up to 1 mas.

3) This GAST is compatible with the IAU 2000 resolutions and must be
   used only in conjunction with other IAU 2000 compatible components
   such as precession-nutation.

4) The result is returned in the range 0 to 2pi.

5) The algorithm is from Capitaine et al. (2003) and IERS Conventions
   2003.

# References

Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
implement the IAU 2000 definition of UT1", Astronomy & Astrophysics,
406, 1135-1149 (2003)

McCarthy, D.D. & Luzum, B.J., "An abridged model of the
precession-nutation of the celestial pole", Celestial Mechanics &
Dynamical Astronomy, 85, 37-49 (2003)

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function gst00b(ut1::AbstractFloat, ut2::AbstractFloat)
    mod2pi(gmst00(ut1, ut2, ut1, ut2) + ee00b(ut1, ut2))
end

"""
    gst06(ut1a::AbstractFloat, ut1b::AbstractFloat, tta::AbstractFloat, ttb::AbstractFloat,
          r::AbstractMatrix{<:AbstractFloat})

Greenwich apparent sidereal time, IAU 2006, given the NPB matrix.

# Input

 - `uta`   -- UT1 as Julian Date (Notes 1,2)
 - `utb`   -- ... as Julian Date
 - `tta`   -- TT as Julian Date (Notes 1,2)
 - `ttb`   -- ... as Julian Date
 - `rnpb`  -- nutation x precession x bias matrix

# Output

 - `gst`   -- Greenwich apparent sidereal time (radians)

# Note

1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
   Julian Dates, apportioned in any convenient way between the
   argument pairs.  For example, JD(UT1)=2450123.7 could be expressed
   in any of these ways, among others:

           uta            utb

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable (in the case of UT; the TT is not at all critical in
   this respect).  The J2000 and MJD methods are good compromises
   between resolution and convenience.  For UT, the date & time method
   is best matched to the algorithm that is used by the Earth rotation
   angle function, called internally: maximum precision is delivered
   when the uta argument is for 0hrs UT1 on the day in question and
   the utb argument lies in the range 0 to 1, or vice versa.

2) Both UT1 and TT are required, UT1 to predict the Earth rotation and
   TT to predict the effects of precession-nutation.  If UT1 is used
   for both purposes, errors of order 100 microarcseconds result.

3) Although the function uses the IAU 2006 series for s+XY/2, it is
   otherwise independent of the precession-nutation model and can in
   practice be used with any equinox-based NPB matrix.

4) The result is returned in the range 0 to 2pi.

# References

Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
"""
function gst06(ut1a::F, ut1b::F, tta::F, ttb::F, r::M) where
   {F<:AbstractFloat, M<:AbstractMatrix{<:AbstractFloat}}
    @inline mod2pi(era00(ut1a, ut1b) - eors(r, s06(tta, ttb, bpn2xy(r)...)))
end

"""
    gst06a(ut1::AbstractFloat, ut2::AbstractFloat, tt1::AbstractFloat, tt2::AbstractFloat)

Greenwich apparent sidereal time (consistent with IAU 2000 and 2006
resolutions).

# Input

 - `uta`   -- UT1 as Julian Date (Notes 1,2)
 - `utb`   -- ... as Julian Date
 - `tta`   -- TT as Julian Date (Notes 1,2)
 - `ttb`   -- ... as Julian Date

# Output

 - `gst`   -- Greenwich apparent sidereal time (radians)

# Note

1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
   Julian Dates, apportioned in any convenient way between the
   argument pairs.  For example, JD(UT1)=2450123.7 could be expressed
   in any of these ways, among others:

           uta            utb

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable (in the case of UT; the TT is not at all critical in
   this respect).  The J2000 and MJD methods are good compromises
   between resolution and convenience.  For UT, the date & time method
   is best matched to the algorithm that is used by the Earth rotation
   angle function, called internally: maximum precision is delivered
   when the uta argument is for 0hrs UT1 on the day in question and
   the utb argument lies in the range 0 to 1, or vice versa.

2) Both UT1 and TT are required, UT1 to predict the Earth rotation and
   TT to predict the effects of precession-nutation.  If UT1 is used
   for both purposes, errors of order 100 microarcseconds result.

3) This GAST is compatible with the IAU 2000/2006 resolutions and must
   be used only in conjunction with IAU 2006 precession and IAU 2000A
   nutation.

4) The result is returned in the range 0 to 2pi.

# References

Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
"""
function gst06a(ut1::AbstractFloat, ut2::AbstractFloat, tt1::AbstractFloat, tt2::AbstractFloat)
    #  Greenwich apparent sidereal time using classical
    #  nutation-precession-bias matrix (IAU 2000A)
    gst06(ut1, ut2, tt1, tt2, pnm06a(tt1, tt2))
end

"""
    gst94(ut1::AbstractFloat, ut2::AbstractFloat)

Greenwich apparent sidereal time (consistent with IAU 1982/94
resolutions).

# Input

 - `uta`   -- UT1 as Julian Date (Notes 1,2)
 - `utb`   -- ... as Julian Date

# Output

 - `gst`   -- Greenwich apparent sidereal time (radians)

# Note

1) The UT1 date uta+utb is a Julian Date, apportioned in any
   convenient way between the argument pair.  For example,
   JD(UT1)=2450123.7 could be expressed in any of these ways, among
   others:

           uta            utb

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 and MJD methods are good compromises between
   resolution and convenience.  For UT, the date & time method is best
   matched to the algorithm that is used by the Earth Rotation Angle
   function, called internally: maximum precision is delivered when
   the uta argument is for 0hrs UT1 on the day in question and the utb
   argument lies in the range 0 to 1, or vice versa.

2) The result is compatible with the IAU 1982 and 1994 resolutions,
   except that accuracy has been compromised for the sake of
   convenience in that UT is used instead of TDB (or TT) to compute
   the equation of the equinoxes.

3) This GAST must be used only in conjunction with contemporaneous IAU
   standards such as 1976 precession, 1980 obliquity and 1982
   nutation.  It is not compatible with the IAU 2000 resolutions.

4) The result is returned in the range 0 to 2pi.

# References

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992)

IAU Resolution C7, Recommendation 3 (1994)
"""
gst94(ut1::AbstractFloat, ut2::AbstractFloat) = mod2pi(gmst82(ut1, ut2) + eqeq94(ut1, ut2))
