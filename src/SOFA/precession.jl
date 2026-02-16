####    Astronomy / Precession-Nutation-Polar Orientation    ####

"""
    bi00()

Frame bias components of IAU 2000 precession-nutation models; part of
the Mathews-Herring-Buffett (MHB2000) nutation series, with additions.

# Output

   dpsibi,depsbi  double  longitude and obliquity corrections
   dra            double  the ICRS RA of the J2000.0 mean equinox

# Note

1) The frame bias corrections in longitude and obliquity (radians) are
   required in order to correct for the offset between the GCRS pole
   and the mean J2000.0 pole.  They define, with respect to the GCRS
   frame, a J2000.0 mean pole that is consistent with the rest of the
   IAU 2000A precession-nutation model.

2) In addition to the displacement of the pole, the complete
   description of the frame bias requires also an offset in right
   ascension.  This is not part of the IAU 2000A model, and is from
   Chapront et al. (2002).  It is returned in radians.

3) This is a supplemented implementation of one aspect of the IAU
   2000A nutation model, formally adopted by the IAU General Assembly
   in 2000, namely MHB2000 (Mathews et al. 2002).

# References

Chapront, J., Chapront-Touze, M. & Francou, G., Astron.  Astrophys.,
387, 700, 2002.

Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation and
precession: New nutation series for nonrigid Earth and insights into
the Earth's interior", J.Geophys.Res., 107, B4, 2002.  The MHB2000
code itself was obtained on 2002 September 9 from
ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
"""
function bi00()
    #  The frame bias corrections in longitude and obliquity, and the ICRS
    #  RA of the J2000.0 equinox (Chapront et al. 2002).
    NamedTuple{(:ψ, :ϵ, :RA)}(
        (deg2rad.(SVector(ψ_bias_2000, ϵ_bias_2000, icrs_ra_2000))/3600.0))
end

"""
    bp00(day1::AbstractFloat, day2::AbstractFloat)

Frame bias and precession, IAU 2000.

# Input

 - `day1`  -- TT as a Julian Date (Note 1)
 - `day2`  -- ... Julian Date (Note 1)

# Output

 - `rb`    -- frame bias matrix (Note 2)
 - `rp`    -- precession matrix (Note 3)
 - `rbp`   -- bias-precession matrix (Note 4)

# Note

1) The TT date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TT)=2450123.7 could be expressed in any of these ways, among
   others:

           date1         date2

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

2) The matrix rb transforms vectors from GCRS to mean J2000.0 by
   applying frame bias.

3) The matrix rp transforms vectors from J2000.0 mean equator and
   equinox to mean equator and equinox of date by applying precession.

4) The matrix rbp transforms vectors from GCRS to mean equator and
   equinox of date by applying frame bias then precession.  It is the
   product rp x rb.

5) It is permissible to re-use the same array in the returned
   arguments.  The arrays are filled in the order given.

# References

"Expressions for the Celestial Intermediate Pole and Celestial
Ephemeris Origin consistent with the IAU 2000A precession- nutation
model", Astron.Astrophys. 400, 1145-1154 (2003)

n.b. The celestial ephemeris origin (CEO) was renamed "celestial
     intermediate origin" (CIO) by IAU 2006 Resolution 2.
"""
function bp00(day1::AbstractFloat, day2::AbstractFloat)
    #  Interval between fundamental epoch J2000.0 and current date (JC).
    Δt = ((day1 - JD2000) + day2)/(100*DAYPERYEAR)
    #  Frame bias
    δψ, δϵ, δra = bi00()
    #  Precession angles (Lieske et al. 1977)
    χA = deg2rad(Polynomial(χ_1977...)(Δt)/3600)
    #  Apply IAU 2000 precession corrections.
    ψA, ωA = deg2rad.((Polynomial(ψ_1977...)(Δt),
                       Polynomial(ω_1977...)(Δt))./3600) .+
                           values(pr00(day1, day2))
    #  Frame bias matrix: GCRS to J2000.0.
    rb = Rx(-δϵ)Ry(δψ*sin(deg2rad(ϵ0_2000/3600)))Rz(δra)
    rp = Rz(χA)Rx(-ωA)Rz(-ψA)Rx(deg2rad(ϵ0_2000/3600))
    #  Bias-precession matrix: GCRS to mean of date.
    (rb = rb, rp = rp, rbp = rp*rb)
end

"""
    bp06(day1::AbstractFloat, day2::AbstractFloat)

Frame bias and precession, IAU 2006.

This function is part of the International Astronomical Union's SOFA
(Standards of Fundamental Astronomy) software collection.

Status:  support function.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... Julian Date

# Output

 - `rb`    -- frame bias matrix (Note 2)
 - `rp`    -- precession matrix (Note 3)
 - `rbp`   -- bias-precession matrix (Note 4)

# Notes

1) The TT date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TT)=2450123.7 could be expressed in any of these ways, among
   others:

           date1         date2

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

2) The matrix rb transforms vectors from GCRS to mean J2000.0 by
   applying frame bias.

3) The matrix rp transforms vectors from mean J2000.0 to mean of date
   by applying precession.

4) The matrix rbp transforms vectors from GCRS to mean of date by
   applying frame bias then precession.  It is the product rp x rb.

5) It is permissible to re-use the same array in the returned
   arguments.  The arrays are filled in the order given.

# References

Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
"""
function bp06(day1::AbstractFloat, day2::AbstractFloat)
    #  B matrix
    rb = fw2m(pfw06(MJD0, MJD00)...)
    #  PxB matrix
    rbp = pmat06(day1, day2)
    #  P matrix
    (rb = rb, rp = rbp*rb', rbp = rbp)
end

"""
    bpn2xy(r::AbstractMatrix{<:AbstractFloat})

Extract from the bias-precession-nutation matrix the X,Y coordinates
of the Celestial Intermediate Pole.

# Input

 - `rbpn`  -- celestial-to-true matrix (Note 1)

# Output

 - `x, y`  -- Celestial Intermediate Pole (Note 2)

# Note

1) The matrix rbpn transforms vectors from GCRS to true equator (and
   CIO or equinox) of date, and therefore the Celestial Intermediate
   Pole unit vector is the bottom row of the matrix.

2) The arguments x,y are components of the Celestial Intermediate Pole
   unit vector in the Geocentric Celestial Reference System.

# References

"Expressions for the Celestial Intermediate Pole and Celestial
Ephemeris Origin consistent with the IAU 2000A precession- nutation
model", Astron.Astrophys. 400, 1145-1154 (2003)

n.b. The celestial ephemeris origin (CEO) was renamed "celestial
     intermediate origin" (CIO) by IAU 2006 Resolution 2.
"""
bpn2xy(r::AbstractMatrix{<:AbstractFloat}) = r[3,1:2]

"""
    c2i00a(day1::AbstractFloat, day2::AbstractFloat)

Form the celestial-to-intermediate matrix for a given date using the
IAU 2000A precession-nutation model.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `rc2i`  -- celestial-to-intermediate matrix (Note 2)

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
   good compromises between resolution and convenience.

2) The matrix rc2i is the first stage in the transformation from
   celestial to terrestrial coordinates:

      [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]

             =  rc2t * [CRS]

   where [CRS] is a vector in the Geocentric Celestial Reference
   System and [TRS] is a vector in the International Terrestrial
   Reference System (see IERS Conventions 2003), ERA is the Earth
   Rotation Angle and RPOM is the polar motion matrix.

3) A faster, but slightly less accurate, result (about 1 mas) can be
   obtained by using instead the eraC2i00b function.

# References

"Expressions for the Celestial Intermediate Pole and Celestial
Ephemeris Origin consistent with the IAU 2000A precession- nutation
model", Astron.Astrophys. 400, 1145-1154 (2003)

n.b. The celestial ephemeris origin (CEO) was renamed "celestial
     intermediate origin" (CIO) by IAU 2006 Resolution 2.

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function c2i00a(day1::F, day2::F) where F<:AbstractFloat
    #  Obtain the celestial-to-true matrix (IAU 2000A) and form the
    #  celestial-to-intermediate matrix
    c2ibpn(day1, day2, pnm00a(day1, day2))
end

"""
    c2i00b(day1::AbstractFloat, day2::AbstractFloat)

Form the celestial-to-intermediate matrix for a given date using the
IAU 2000B precession-nutation model.

# Input

 - `day1`  -- TT as 2-part Julian Date (Note 1)
 - `day2`  -- ... Julian Date

# Output

 - `rc2i`  -- celestial-to-intermediate matrix (Note 2)

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
   good compromises between resolution and convenience.

2) The matrix rc2i is the first stage in the transformation from
   celestial to terrestrial coordinates:

      [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]

             =  rc2t * [CRS]

   where [CRS] is a vector in the Geocentric Celestial Reference
   System and [TRS] is a vector in the International Terrestrial
   Reference System (see IERS Conventions 2003), ERA is the Earth
   Rotation Angle and RPOM is the polar motion matrix.

3) The present function is faster, but slightly less accurate (about 1
   mas), than the eraC2i00a function.

# References

"Expressions for the Celestial Intermediate Pole and Celestial
Ephemeris Origin consistent with the IAU 2000A precession- nutation
model", Astron.Astrophys. 400, 1145-1154 (2003)

n.b. The celestial ephemeris origin (CEO) was renamed "celestial
     intermediate origin" (CIO) by IAU 2006 Resolution 2.

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function c2i00b(day1::F, day2::F) where F<:AbstractFloat
    #  Obtain the celestial-to-true matrix (IAU 2000B) and form the
    #  celestial-to-intermediate matrix
    c2ibpn(day1, day2, pnm00b(day1, day2))
end

"""
    c2i06a(day1::AbstractFloat, day2::AbstractFloat)

Form the celestial-to-intermediate matrix for a given date using the
IAU 2006 precession and IAU 2000A nutation models.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... Julian Date

# Output

 - `rc2i`  -- celestial-to-intermediate matrix (Note 2)

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
   good compromises between resolution and convenience.

2) The matrix rc2i is the first stage in the transformation from
   celestial to terrestrial coordinates:

      [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]

             =  RC2T * [CRS]

   where [CRS] is a vector in the Geocentric Celestial Reference
   System and [TRS] is a vector in the International Terrestrial
   Reference System (see IERS Conventions 2003), ERA is the Earth
   Rotation Angle and RPOM is the polar motion matrix.

# References

McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003), IERS
Technical Note No. 32, BKG
"""
function c2i06a(day1::F, day2::F) where F<:AbstractFloat
    #  Obtain the celestial-to-true matrix (IAU 2006/2000A), extract
    #  x, y coordinates.
    x, y = bpn2xy(pnm06a(day1, day2))
    #  Obtain the CIO locator, and form the celestial-to-intermediate matrix
    c2ixys(x, y, s06(day1, day2, x, y))
end

"""
    c2ibpn(day1::AbstractFloat, day2::AbstractFloat, r::AbstractMatrix{<:AbstractFloat})

Form the celestial-to-intermediate matrix for a given date given the
bias-precession-nutation matrix.  IAU 2000.

#Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... Julian Date
 - `rbpn`  -- celestial-to-true matrix (Note 2)

# Output

 - `rc2i`  -- celestial-to-intermediate matrix (Note 3)

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
   good compromises between resolution and convenience.

2) The matrix rbpn transforms vectors from GCRS to true equator (and
   CIO or equinox) of date.  Only the CIP (bottom row) is used.

3) The matrix rc2i is the first stage in the transformation from
   celestial to terrestrial coordinates:

      [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]

            = RC2T * [CRS]

   where [CRS] is a vector in the Geocentric Celestial Reference
   System and [TRS] is a vector in the International Terrestrial
   Reference System (see IERS Conventions 2003), ERA is the Earth
   Rotation Angle and RPOM is the polar motion matrix.

4) Although its name does not include "00", This function is in fact
   specific to the IAU 2000 models.

# References

"Expressions for the Celestial Intermediate Pole and Celestial
Ephemeris Origin consistent with the IAU 2000A precession- nutation
model", Astron.Astrophys. 400, 1145-1154 (2003)

n.b. The celestial ephemeris origin (CEO) was renamed "celestial
     intermediate origin" (CIO) by IAU 2006 Resolution 2.

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function c2ibpn(day1::F, day2::F, r::M) where
   {F<:AbstractFloat, M<:AbstractMatrix{<:AbstractFloat}}
    c2ixy(day1, day2, bpn2xy(r)...)
end

"""
    c2ixy(day1::AbstractFloat, day2::AbstractFloat, x::AbstractFloat, y::AbstractFloat)

Form the celestial to intermediate-frame-of-date matrix for a given
date when the CIP X,Y coordinates are known.  IAU 2000.

#Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... Julian Date
 - `x, y`  -- Celestial Intermediate Pole (Note 2)

# Output

 - `rc2i`  -- celestial-to-intermediate matrix (Note 3)

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
   good compromises between resolution and convenience.

2) The Celestial Intermediate Pole coordinates are the x,y components
   of the unit vector in the Geocentric Celestial Reference System.

3) The matrix rc2i is the first stage in the transformation from
   celestial to terrestrial coordinates:

      [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]

            = RC2T * [CRS]

   where [CRS] is a vector in the Geocentric Celestial Reference
   System and [TRS] is a vector in the International Terrestrial
   Reference System (see IERS Conventions 2003), ERA is the Earth
   Rotation Angle and RPOM is the polar motion matrix.

4) Although its name does not include "00", This function is in fact
   specific to the IAU 2000 models.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function c2ixy(day1::F, day2::F, x::F, y::F) where F<:AbstractFloat
    c2ixys(x, y, s00(day1, day2, x, y))
end

"""
    c2ixys(x::AbstractFloat, y::AbstractFloat, s::AbstractFloat)

Form the celestial to intermediate-frame-of-date matrix given the CIP
X,Y and the CIO locator s.

# Input

 - `x, y`  -- Celestial Intermediate Pole (Note 1)
 - `s`     -- The CIO locator s (Note 2)

# Output

 - `rc2i`  -- celestial-to-intermediate matrix (Note 3)

# Note

1) The Celestial Intermediate Pole coordinates are the x,y components
   of the unit vector in the Geocentric Celestial Reference System.

2) The CIO locator s (in radians) positions the Celestial Intermediate
   Origin on the equator of the CIP.

3) The matrix rc2i is the first stage in the transformation from
   celestial to terrestrial coordinates:

      [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]

            = RC2T * [CRS]

   where [CRS] is a vector in the Geocentric Celestial Reference
   System and [TRS] is a vector in the International Terrestrial
   Reference System (see IERS Conventions 2003), ERA is the Earth
   Rotation Angle and RPOM is the polar motion matrix.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function c2ixys(x::F, y::F, s::F) where F<:AbstractFloat
    r = x*x + y*y
    e = r > 0.0 ? atan(y, x) : 0.0
    Rz(-(e+s))Ry(atan(sqrt(r/(1.0 - r))))Rz(e)
end

"""
    c2t00a(tt1::AbstractFloat, tt2::AbstractFloat, ut1::AbstractFloat, ut2::AbstractFloat, x::AbstractFloat,
           y::AbstractFloat)

Form the celestial to terrestrial matrix given the date, the UT1 and
the polar motion, using the IAU 2000A precession-nutation model.

# Input

 - `tta, ttb` -- TT as a 2-part Julian Date (Note 1)
 - `uta, utb` -- UT1 as a 2-part Julian Date (Note 1)
 - `xp, yp`   -- CIP coordinates (radians, Note 2)

# Outpu

 - `rc2t`  -- celestial-to-terrestrial matrix (Note 3)

# Note

1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
   apportioned in any convenient way between the arguments uta and
   utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
   these ways, among others:

           uta            utb

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 and MJD methods are good compromises between
   resolution and convenience.  In the case of uta,utb, the date &
   time method is best matched to the Earth rotation angle algorithm
   used: maximum precision is delivered when the uta argument is for
   0hrs UT1 on the day in question and the utb argument lies in the
   range 0 to 1, or vice versa.

2) The arguments xp and yp are the coordinates (in radians) of the
   Celestial Intermediate Pole with respect to the International
   Terrestrial Reference System (see IERS Conventions 2003), measured
   along the meridians 0 and 90 deg west respectively.

3) The matrix rc2t transforms from celestial to terrestrial
   coordinates:

      [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]

            = rc2t * [CRS]

   where [CRS] is a vector in the Geocentric Celestial Reference
   System and [TRS] is a vector in the International Terrestrial
   Reference System (see IERS Conventions 2003), RC2I is the
   celestial-to-intermediate matrix, ERA is the Earth rotation angle
   and RPOM is the polar motion matrix.

4) A faster, but slightly less accurate, result (about 1 mas) can be
   obtained by using instead the eraC2t00b function.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function c2t00a(tt1::F, tt2::F, ut1::F, ut2::F, x::F, y::F) where F<:AbstractFloat
    #  Form the celestial-to-intermediate matrix for this TT (IAU 2000A),
    #  predict the Earth rotation angle for this UT1, estimate s', form the
    #  polar motion matrix, and combine to form the celestial-to-terrestrial
    #  matrix.
    c2tcio(c2i00a(tt1, tt2), era00(ut1, ut2), pom00(x, y, sp00(tt1, tt2)))
end

"""
    c2t00b(tt1::AbstractFloat, tt2::AbstractFloat, ut1::AbstractFloat, ut2::AbstractFloat, x::AbstractFloat,
           y::AbstractFloat)

Form the celestial to terrestrial matrix given the date, the UT1 and
the polar motion, using the IAU 2000B precession-nutation model.

# Input

 - `tta, ttb` -- TT as a 2-part Julian Date (Note 1)
 - `uta, utb` -- UT1 as a 2-part Julian Date (Note 1)
 - `xp, yp`   -- coordinates of the pole (radians, Note 2)

# Output

 - `rc2t`  -- celestial-to-terrestrial matrix (Note 3)

# Note

1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
   apportioned in any convenient way between the arguments uta and
   utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
   these ways, among others:

           uta            utb

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 and MJD methods are good compromises between
   resolution and convenience.  In the case of uta,utb, the date &
   time method is best matched to the Earth rotation angle algorithm
   used: maximum precision is delivered when the uta argument is for
   0hrs UT1 on the day in question and the utb argument lies in the
   range 0 to 1, or vice versa.

2) The arguments xp and yp are the coordinates (in radians) of the
   Celestial Intermediate Pole with respect to the International
   Terrestrial Reference System (see IERS Conventions 2003), measured
   along the meridians 0 and 90 deg west respectively.

3) The matrix rc2t transforms from celestial to terrestrial
   coordinates:

      [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]

            = rc2t * [CRS]

   where [CRS] is a vector in the Geocentric Celestial Reference
   System and [TRS] is a vector in the International Terrestrial
   Reference System (see IERS Conventions 2003), RC2I is the
   celestial-to-intermediate matrix, ERA is the Earth rotation angle
   and RPOM is the polar motion matrix.

4) The present function is faster, but slightly less accurate (about 1
   mas), than the eraC2t00a function.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function c2t00b(tt1::F, tt2::F, ut1::F, ut2::F, x::F, y::F) where F<:AbstractFloat
    #  Form the celestial-to-intermediate matrix for this TT (IAU 2000A),
    #  predict the Earth rotation angle for this UT1, estimate s', form the
    #  polar motion matrix, and combine to form the celestial-to-terrestrial
    #  matrix.
    c2tcio(c2i00b(tt1, tt2), era00(ut1, ut2), pom00(x, y, 0.0))
end

"""
    c2t06a(tt1::AbstractFloat, tt2::AbstractFloat, ut1::AbstractFloat, ut2::AbstractFloat, x::AbstractFloat,
           y::AbstractFloat)

Form the celestial to terrestrial matrix given the date, the UT1 and
the polar motion, using the IAU 2006/2000A precession-nutation model.

#  Input

 - `tta, ttb` -- TT as a 2-part Julian Date (Note 1)
 - `uta, utb` -- UT1 as a 2-part Julian Date (Note 1)
 - `xp, yp`   -- coordinates of the pole (radians, Note 2)

# Output

 - `rc2t`   -- celestial-to-terrestrial matrix (Note 3)

# Note

1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
   apportioned in any convenient way between the two arguments uta and
   utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
   these ways, among others:

           uta            utb

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 and MJD methods are good compromises between
   resolution and convenience.  In the case of uta,utb, the date &
   time method is best matched to the Earth rotation angle algorithm
   used: maximum precision is delivered when the uta argument is for
   0hrs UT1 on the day in question and the utb argument lies in the
   range 0 to 1, or vice versa.

   where [CRS] is a vector in the Geocentric Celestial Reference
   System and [TRS] is a vector in the International Terrestrial
   Reference System (see IERS Conventions 2003), RC2I is the
   celestial-to-intermediate matrix, ERA is the Earth rotation angle
   and RPOM is the polar motion matrix.

# References

McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003), IERS
Technical Note No. 32, BKG
"""
function c2t06a(tt1::F, tt2::F, ut1::F, ut2::F, x::F, y::F) where F<:AbstractFloat
    #  Form the celestial-to-intermediate matrix for this TT (IAU 2000A),
    #  predict the Earth rotation angle for this UT1, estimate s', form the
    #  polar motion matrix, and combine to form the celestial-to-terrestrial
    #  matrix.
    c2tcio(c2i06a(tt1, tt2), era00(ut1, ut2), pom00(x, y, sp00(tt1, tt2)))
end

"""
    c2tcio(c2i::AbstractMatrix{<:AbstractFloat}, era::AbstractFloat, pm::AbstractMatrix{<:AbstractFloat}) = pm*Rz(era)*c2i

Assemble the celestial to terrestrial matrix from CIO-based components
(the celestial-to-intermediate matrix, the Earth Rotation Angle and
the polar motion matrix).

# Input

 - `rc2i`  -- celestial-to-intermediate matrix
 - `era`   -- Earth rotation angle (radians)
 - `rpom`  -- polar-motion matrix

# Output

 - `rc2t`  -- celestial-to-terrestrial matrix

# Note

1) This function constructs the rotation matrix that transforms
   vectors in the celestial system into vectors in the terrestrial
   system.  It does so starting from precomputed components, namely
   the matrix which rotates from celestial coordinates to the
   intermediate frame, the Earth rotation angle and the polar motion
   matrix.  One use of the present function is when generating a
   series of celestial-to-terrestrial matrices where only the Earth
   Rotation Angle changes, avoiding the considerable overhead of
   recomputing the precession-nutation more often than necessary to
   achieve given accuracy objectives.

2) The relationship between the arguments is as follows:

      [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]

            = rc2t * [CRS]

   where [CRS] is a vector in the Geocentric Celestial Reference
   System and [TRS] is a vector in the International Terrestrial
   Reference System (see IERS Conventions 2003).

# References

McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003), IERS
Technical Note No. 32, BKG
"""
function c2tcio(c2i::M, era::F, pm::M) where
   {M<:AbstractMatrix{<:AbstractFloat}, F<:AbstractFloat}
   pm*Rz(era)*c2i
end

"""
    c2teqx(bpn::AbstractMatrix{<:AbstractFloat}, gst::AbstractFloat, pm::AbstractMatrix{<:AbstractFloat}) = pm*Rz(gst)*bpn

Assemble the celestial to terrestrial matrix from equinox-based
components (the celestial-to-true matrix, the Greenwich Apparent
Sidereal Time and the polar motion matrix).

# Input

 - `rbpn`  -- celestial-to-true matrix
 - `gst`   -- Greenwich (apparent) Sidereal Time (radians)
 - `rpom`  -- polar-motion matrix

# Output

 - `rc2t`  -- celestial-to-terrestrial matrix (Note 2)

# Note

1) This function constructs the rotation matrix that transforms
   vectors in the celestial system into vectors in the terrestrial
   system.  It does so starting from precomputed components, namely
   the matrix which rotates from celestial coordinates to the true
   equator and equinox of date, the Greenwich Apparent Sidereal Time
   and the polar motion matrix.  One use of the present function is
   when generating a series of celestial-to-terrestrial matrices where
   only the Sidereal Time changes, avoiding the considerable overhead
   of recomputing the precession-nutation more often than necessary to
   achieve given accuracy objectives.

2) The relationship between the arguments is as follows:

      [TRS] = rpom * R_3(gst) * rbpn * [CRS]

            = rc2t * [CRS]

   where [CRS] is a vector in the Geocentric Celestial Reference
   System and [TRS] is a vector in the International Terrestrial
   Reference System (see IERS Conventions 2003).

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function c2teqx(bpn::M, gst::F, pm::M) where
   {M<:AbstractMatrix{<:AbstractFloat}, F<:AbstractFloat}
   pm*Rz(gst)*bpn
end

"""
    c2tpe(tt1::AbstractFloat, tt2::AbstractFloat, ut1::AbstractFloat, ut2::AbstractFloat, ψ::AbstractFloat,
          ϵ::AbstractFloat, xp::AbstractFloat, yp::AbstractFloat)

Form the celestial to terrestrial matrix given the date, the UT1, the
nutation and the polar motion.  IAU 2000.

# Input

 - `tta, ttb`   -- TT as a 2-part Julian Date (Note 1)
 - `uta, utb`   -- UT1 as a 2-part Julian Date (Note 1)
 - `dpsi, deps` -- nutation (Note 2)
 - `xp, yp`     -- coordinates of the pole (radians, Note 3)

# Output

 - `rc2t`  -- celestial-to-terrestrial matrix (Note 4)

# Note

1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
   apportioned in any convenient way between the arguments uta and
   utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
   these ways, among others:

           uta            utb

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 and MJD methods are good compromises between
   resolution and convenience.  In the case of uta,utb, the date &
   time method is best matched to the Earth rotation angle algorithm
   used: maximum precision is delivered when the uta argument is for
   0hrs UT1 on the day in question and the utb argument lies in the
   range 0 to 1, or vice versa.

2) The caller is responsible for providing the nutation components;
   they are in longitude and obliquity, in radians and are with
   respect to the equinox and ecliptic of date.  For high-accuracy
   applications, free core nutation should be included as well as any
   other relevant corrections to the position of the CIP.

3) The arguments xp and yp are the coordinates (in radians) of the
   Celestial Intermediate Pole with respect to the International
   Terrestrial Reference System (see IERS Conventions 2003), measured
   along the meridians 0 and 90 deg west respectively.

4) The matrix rc2t transforms from celestial to terrestrial
   coordinates:

      [TRS] = RPOM * R_3(GST) * RBPN * [CRS]

            = rc2t * [CRS]

   where [CRS] is a vector in the Geocentric Celestial Reference
   System and [TRS] is a vector in the International Terrestrial
   Reference System (see IERS Conventions 2003), RBPN is the
   bias-precession-nutation matrix, GST is the Greenwich (apparent)
   Sidereal Time and RPOM is the polar motion matrix.

5) Although its name does not include "00", This function is in fact
   specific to the IAU 2000 models.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function c2tpe(tt1::F, tt2::F, ut1::F, ut2::F, ψ::F, ϵ::F, xp::F, yp::F) where
      F<:AbstractFloat
    #  Form the celestial-to-intermediate matrix for this TT
    ϵA, rb, rp, rbp, rn, rbpn = values(pn00(tt1, tt2, ψ, ϵ))
    #  Predict the Greenwich Mean Sidereal Time for this UT1 and TT, predict
    #  the equation of the equinoxes for this TT and nutation, estimate s',
    #  form the polar motion matrix, and combine to form the
    #  celestial-to-terrestrial matrix.
    c2teqx(rbpn, gmst00(ut1, ut2, tt1, tt2) + ee00(tt1, tt2, ϵA, ψ),
           pom00(xp, yp, sp00(tt1, tt2)))
end

"""
    c2txy(tt1::AbstractFloat, tt2::AbstractFloat, ut1::AbstractFloat, ut2::AbstractFloat, x, y, xp, yp)

Form the celestial to terrestrial matrix given the date, the UT1, the
CIP coordinates and the polar motion.  IAU 2000.

# Output

 - `tta, ttb` -- TT as a 2-part Julian Date (Note 1)
 - `uta, utb` -- UT1 as a 2-part Julian Date (Note 1)
 - `x, y`     -- Celestial Intermediate Pole (Note 2)
 - `xp, yp`   -- coordinates of the pole (radians, Note 3)

# Output

 - `rc2t`  -- celestial-to-terrestrial matrix (Note 4)

# Note

1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
   apportioned in any convenient way between the arguments uta and
   utb.  For example, JD(UT1)=2450123.7 could be expressed in any o
   these ways, among others:

           uta            utb

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 and MJD methods are good compromises between
   resolution and convenience.  In the case of uta,utb, the date &
   time method is best matched to the Earth rotation angle algorithm
   used: maximum precision is delivered when the uta argument is for
   0hrs UT1 on the day in question and the utb argument lies in the
   range 0 to 1, or vice versa.

2) The Celestial Intermediate Pole coordinates are the x,y components
   of the unit vector in the Geocentric Celestial Reference System.

3) The arguments xp and yp are the coordinates (in radians) of the
   Celestial Intermediate Pole with respect to the International
   Terrestrial Reference System (see IERS Conventions 2003), measured
   along the meridians 0 and 90 deg west respectively.

4) The matrix rc2t transforms from celestial to terrestrial
   coordinates:

      [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]

            = rc2t * [CRS]

   where [CRS] is a vector in the Geocentric Celestial Reference
   System and [TRS] is a vector in the International Terrestrial
   Reference System (see IERS Conventions 2003), ERA is the Earth
   Rotation Angle and RPOM is the polar motion matrix.

5) Although its name does not include "00", This function is in fact
   specific to the IAU 2000 models.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function c2txy(tt1::F, tt2::F, ut1::F, ut2::F, x, y, xp, yp) where F<:AbstractFloat
    #  Form the celestial-to-intermediate matrix for this TT, predict the Earth
    #  rotation angle for this UT1, estimate s', form the polar motion matrix, and
    #  combine to form the celestial-to-terrestrial matrix.
    c2tcio(c2ixy(tt1, tt2, x, y), era00(ut1, ut2), pom00(xp, yp, sp00(tt1, tt2)))
end

"""
    eo06a(day1::AbstractFloat, day2::AbstractFloat)

Equation of the origins, IAU 2006 precession and IAU 2000A nutation.

# Output

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `eos`   -- the equation of the origins in radians

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
   good compromises between resolution and convenience.

2) The equation of the origins is the distance between the true
   equinox and the celestial intermediate origin and, equivalently,
   the difference between Earth rotation angle and Greenwich apparent
   sidereal time (ERA-GST).  It comprises the precession (since
   J2000.0) in right ascension plus the equation of the equinoxes
   (including the small correction terms).

# References

Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
"""
function eo06a(day1::F, day2::F) where F<:AbstractFloat
    #  Classical nutation-precession-bias matrix.
    bpn = pnm06a(day1, day2)
    #  Extract CIP coordinates, the CIO locator (s), and solve for
    #  the equation of the origin.
    eors(bpn, s06(day1, day2, bpn2xy(bpn)...))
end

"""
    eors(r::AbstractMatrix{<:AbstractFloat}, s::AbstractFloat)

Equation of the origins, given the classical NPB matrix and the
quantity s.

# Input

 - `rnpb`  -- classical nutation x precession x bias matrix
 - `s`     -- the quantity s (the CIO locator) in radians

# Output

 - `eos`   -- the equation of the origins in radians

# Note

1) The equation of the origins is the distance between the true
    equinox and the celestial intermediate origin and, equivalently,
    the difference between Earth rotation angle and Greenwich apparent
    sidereal time (ERA-GST).  It comprises the precession (since
    J2000.0) in right ascension plus the equation of the equinoxes
    (including the small correction terms).

2)  The algorithm is from Wallace & Capitaine (2006).

# References

Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

Wallace, P. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
"""
function eors(r::M, s::F) where
   {M<:AbstractMatrix{<:AbstractFloat}, F<:AbstractFloat}
    #  Evaluate Wallace & Capitaine (2006) expression (16).
    v = r*SVector(1.0 - r[3,1]^2/(1.0 + r[3,3]),
      -r[3,2]*r[3,1]/(1.0 + r[3,3]), -r[3,1])
    v[1] != 0 || v[2] != 0 ? s - atan(v[2], v[1]) : s
end

"""
    fw2m(γ::AbstractFloat, ϕ::AbstractFloat, ψ::AbstractFloat, ϵ::AbstractFloat)

Form rotation matrix given the Fukushima-Williams angles.

# Input

 - `γb`    -- F-W angle γ_bar (radians)
 - `ϕb`    -- F-W angle ϕ_bar (radians)
 - `ψ`     -- F-W angle ψ (radians)
 - `ϵ`     -- F-W angle ϵ (radians)

# Output

 - `r`     -- rotation matrix

# Note

1) Naming the following points:

         e = J2000.0 ecliptic pole,
         p = GCRS pole,
         E = ecliptic pole of date,
   and   P = CIP,

   the four Fukushima-Williams angles are as follows:

      γb = γ = epE
      ϕb = ϕ = pE
      ψ  = ψ = pEP
      ϵ  = ϵ = EP

2) The matrix representing the combined effects of frame bias,
   precession and nutation is:

      NxPxB = R_1(-ϵ).R_3(-ψ).R_1(ϕb).R_3(γb)

3) The present function can construct three different matrices,
   depending on which angles are supplied as the arguments gamb, ϕb, ψ
   and ϵ:

   o To obtain the nutation x precession x frame bias matrix, first
      generate the four precession angles known conventionally as
      γ_bar, ϕ_bar, ψ_bar and ϵ_A, then generate the nutation
      components Dψ and Dϵ and add them to ψ_bar and ϵ_A, and finally
      call the present function using those four angles as arguments.

   o To obtain the precession x frame bias matrix, generate the four
      precession angles and call the present function.

   o To obtain the frame bias matrix, generate the four precession
      angles for date J2000.0 and call the present function. The
      nutation-only and precession-only matrices can if necessary be
      obtained by combining these three appropriately.

# References

Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
"""
function fw2m(γ::F, ϕ::F, ψ::F, ϵ::F)  where F<:AbstractFloat
   Rx(-ϵ)Rz(-ψ)Rx(ϕ)Rz(γ)
end

"""
    fw2xy(γ::AbstractFloat, ϕ::AbstractFloat, ψ::AbstractFloat, ϵ::AbstractFloat)

CIP X,Y given Fukushima-Williams bias-precession-nutation angles.

# Input

 - `γb`    -- F-W angle gamma_bar (radians)
 - `ϕb`    -- F-W angle phi_bar (radians)
 - `ψ`     -- F-W angle psi (radians)
 - `ϵ`     -- F-W angle epsilon (radians)

# Output

 - `x, y`  -- CIP unit vector X,Y

# Note

1) Naming the following points:

         e = J2000.0 ecliptic pole,
         p = GCRS pole
         E = ecliptic pole of date,
   and   P = CIP,

   the four Fukushima-Williams angles are as follows:

      γb = γ = epE
      ϕb = ϕ = pE
      ψ = ψ = pEP
      ϵ = ϵ = EP

2) The matrix representing the combined effects of frame bias,
   precession and nutation is:

      NxPxB = R_1(-ϵA).R_3(-ψ).R_1(ϕb).R_3(γb)

   The returned values x,y are elements [2][0] and [2][1] of the
   matrix.  Near J2000.0, they are essentially angles in radians.

# References

Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
"""
function fw2xy(γ::AbstractFloat, ϕ::AbstractFloat, ψ::AbstractFloat, ϵ::AbstractFloat)
    (Rx(-ϵ)Rz(-ψ)Rx(ϕ)Rz(γ))[3,1:2]
end

"""
    ltp(epoch::AbstractFloat)

Long-term precession matrix.

# Input

- `epoch` -- Julian epoch (TT)

# Output

- `r`     -- precession matrix, J2000.0 to date

# Note

1) The matrix is in the sense

      P_date = rp x P_J2000,

   where P_J2000 is a vector with respect to the J2000.0 mean equator
   and equinox and P_date is the same vector with respect to the
   equator and equinox of epoch epj.

2) The Vondrak et al. (2011, 2012) 400 millennia precession model
   agrees with the IAU 2006 precession at J2000.0 and stays within 100
   microarcseconds during the 20th and 21st centuries.  It is accurate
   to a few arcseconds throughout the historical period, worsening to
   a few tenths of a degree at the end of the +/- 200,000 year time
   span.

# References

Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
expressions, valid for long time intervals, Astron.Astrophys. 534, A22

Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
expressions, valid for long time intervals (Corrigendum),
Astron.Astrophys. 541, C1
"""
function ltp(epoch::AbstractFloat)
    #  Equatorial and ecliptic poles
    @inline equ, ecl = ltpequ(epoch), ltpecl(epoch)

    #  Create matrix
    eqx = vec2mat(equ)*ecl/norm(vec2mat(equ)*ecl)
    vcat(eqx', (vec2mat(equ)*eqx)', equ')
end

"""
    ltpb(epoch::AbstractFloat)

Long-term precession matrix, including ICRS frame bias.

# Input

 - `epoch` -- Julian epoch (TT)

# Output

 - `r`     -- precession-bias matrix, J2000.0 to date

# Note

1) The matrix is in the sense

      P_date = rpb x P_ICRS,

   where P_ICRS is a vector in the Geocentric Celestial Reference
   System, and P_date is the vector with respect to the Celestial
   Intermediate Reference System at that date but with nutation
   neglected.

2) A first order frame bias formulation is used, of sub-
   microarcsecond accuracy compared with a full 3D rotation.

3) The Vondrak et al. (2011, 2012) 400 millennia precession model
   agrees with the IAU 2006 precession at J2000.0 and stays within 100
   microarcseconds during the 20th and 21st centuries.  It is accurate
   to a few arcseconds throughout the historical period, worsening to
   a few tenths of a degree at the end of the +/- 200,000 year time
   span.

# References

Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
expressions, valid for long time intervals, Astron.Astrophys. 534, A22

Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
expressions, valid for long time intervals (Corrigendum),
Astron.Astrophys. 541, C1
"""
function ltpb(epoch::AbstractFloat)
    #  Apply frame bias
    @inline ltp(epoch)*(
        SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0) .+ deg2rad(1/3600)*
        SMatrix{3,3}(0.0, -dα0_2010, ϵ0_2010, dα0_2010, 0.0, η0_2010,
         -ϵ0_2010, -η0_2010, 0.0))
end

"""
    ltpecl(epoch::AbstractFloat)

Long-term precession of the ecliptic.

# Input

 - `epoch` -- Julian epoch (TT)

# Output

 - `pos`   -- ecliptic pole unit vector

# Note

1) The returned vector is with respect to the J2000.0 mean equator and
   equinox.

2) The Vondrak et al. (2011, 2012) 400 millennia precession model
   agrees with the IAU 2006 precession at J2000.0 and stays within 100
   microarcseconds during the 20th and 21st centuries.  It is accurate
   to a few arcseconds throughout the historical period, worsening to
   a few tenths of a degree at the end of the +/- 200,000 year time
   span.

# References

Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
expressions, valid for long time intervals, Astron.Astrophys. 534, A22

Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
expressions, valid for long time intervals (Corrigendum),
Astron.Astrophys. 541, C1
"""
function ltpecl(epoch::AbstractFloat)
    #  Centuries since J2000
    Δt = (epoch - 2000.0)/100.0
    ϕ  = 2π*Δt./ecl_ϕ_2011
    p, q = deg2rad(1/3600)*SVector(
        Polynomial(ecl_pA_0_2011...)(Δt) +
        sum(cos.(ϕ).*ecl_pA_c_2011 .+ sin.(ϕ).*ecl_pA_s_2011),
        Polynomial(ecl_qA_0_2011...)(Δt) +
        sum(cos.(ϕ).*ecl_qA_c_2011 .+ sin.(ϕ).*ecl_qA_s_2011))
    one_pp_qq = 1 - p*p - q*q
    w, ϵ0 = one_pp_qq < 0.0 ? 0.0 : sqrt(one_pp_qq), deg2rad(1/3600)*ϵ0_2006
    s, c = sincos(ϵ0)
    return SVector(p, -(q * c + w * s), -(q * s - w * c))
end

"""
    ltpequ(epoch::AbstractFloat)

Long-term precession of the equator.

# Input

 - `epoch` -- Julian epoch (TT)

# Output

 - `pos`   -- equator pole unit vector

# Note

1) The returned vector is with respect to the J2000.0 mean equator and
   equinox.

2) The Vondrak et al. (2011, 2012) 400 millennia precession model
   agrees with the IAU 2006 precession at J2000.0 and stays within 100
   microarcseconds during the 20th and 21st centuries.  It is accurate
   to a few arcseconds throughout the historical period, worsening to
   a few tenths of a degree at the end of the +/- 200,000 year time
   span.

# References

Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
expressions, valid for long time intervals, Astron.Astrophys. 534, A22

Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
expressions, valid for long time intervals (Corrigendum),
Astron.Astrophys. 541, C1
"""
function ltpequ(epoch::AbstractFloat)
    #  Centuries since J2000
    Δt = (epoch - 2000.0)/100.0
    ϕ  = 2π*Δt./equ_ϕ_2011
    x, y = deg2rad(1/3600)*SVector(
        Polynomial(equ_xA_0_2011...)(Δt) +
        sum(cos.(ϕ).*equ_xA_c_2011 .+ sin.(ϕ).*equ_xA_s_2011),
        Polynomial(equ_yA_0_2011...)(Δt) +
        sum(cos.(ϕ).*equ_yA_c_2011 .+ sin.(ϕ).*equ_yA_s_2011))
    SVector(x, y, (1 - x*x - y*y) < 0.0 ? 0.0 : sqrt(1 - x*x - y*y))
end

"""
    num00a(day1::AbstractFloat, day2::AbstractFloat)

Form the matrix of nutation for a given date, IAU 2000A model.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `rmatn` -- nutation matrix

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
   good compromises between resolution and convenience.

2) The matrix operates in the sense V(true) = rmatn * V(mean), where
   the p-vector V(true) is with respect to the true equatorial triad
   of date and the p-vector V(mean) is with respect to the mean
   equatorial triad of date.

3) A faster, but slightly less accurate, result (about 1 mas) can be
   obtained by using instead the eraNum00b function.

# References

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992), Section 3.222-3
(p114).
"""
num00a(day1::AbstractFloat, day2::AbstractFloat) = pn00a(day1, day2)[:rn]

"""
    num00b(day1::AbstractFloat, day2::AbstractFloat)

Form the matrix of nutation for a given date, IAU 2000B model.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... Julian Date

# Output

 - `rmatn` -- nutation matrix

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
   good compromises between resolution and convenience.

2) The matrix operates in the sense V(true) = rmatn * V(mean), where
   the p-vector V(true) is with respect to the true equatorial triad
   of date and the p-vector V(mean) is with respect to the mean
   equatorial triad of date.

3) The present function is faster, but slightly less accurate (about 1
   mas), than the eraNum00a function.

# References

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992), Section 3.222-3
(p114).
"""
num00b(day1::AbstractFloat, day2::AbstractFloat) = pn00b(day1, day2)[:rn]

"""
    num06a(day1::AbstractFloat, day2::AbstractFloat)

Form the matrix of nutation for a given date, IAU 2006/2000A model.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... Julian Date

# Output

 - `rmatn` -- nutation matrix

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
   good compromises between resolution and convenience.

2) The matrix operates in the sense V(true) = rmatn * V(mean), where
   the p-vector V(true) is with respect to the true equatorial triad
   of date and the p-vector V(mean) is with respect to the mean
   equatorial triad of date.
"""
function num06a(day1::AbstractFloat, day2::AbstractFloat)
    #  Mean obliquity, nutation components, and nutation matrix
    numat(obl06(day1, day2), nut06a(day1, day2)...)
end

"""
    numat(ϵA::AbstractFloat, δψ::AbstractFloat, δϵ::AbstractFloat)

Form the matrix of nutation.

# Input

 - `epsa`  -- mean obliquity of date (Note 1)
 - `dpsi, deps` -- nutation (Note 2)

# Output

 - `rmatn` -- nutation matrix (Note 3)

# Note

1) The supplied mean obliquity epsa, must be consistent with the
   precession-nutation models from which dpsi and deps were obtained.

2) The caller is responsible for providing the nutation components;
   they are in longitude and obliquity, in radians and are with
   respect to the equinox and ecliptic of date.

3) The matrix operates in the sense V(true) = rmatn * V(mean), where
   the p-vector V(true) is with respect to the true equatorial triad
   of date and the p-vector V(mean) is with respect to the mean
   equatorial triad of date.

# References

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992), Section 3.222-3
(p114).
"""
function numat(ϵA::AbstractFloat, δψ::AbstractFloat, δϵ::AbstractFloat)
    Rx(-(ϵA + δϵ))Rz(-δψ)Rx(ϵA)
end

"""
    nut00a(day1::AbstractFloat, day2::AbstractFloat)

Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation
with free core nutation omitted).

# Input

 - `day1`   -- TT as Julian Date (Note 1)
 - `day2`   -- ... as Julian Date

# Ouput

 - `dpsi, deps` -- nutation, luni-solar + planetary (Note 2)

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
   good compromises between resolution and convenience.

2) The nutation components in longitude and obliquity are in radians
   and with respect to the equinox and ecliptic of date.  The
   obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
   value of 84381.448 arcsec.

   Both the luni-solar and planetary nutations are included.  The
   latter are due to direct planetary nutations and the perturbations
   of the lunar and terrestrial orbits.

3) The function computes the MHB2000 nutation series with the
   associated corrections for planetary nutations.  It is an
   implementation of the nutation part of the IAU 2000A precession-
   nutation model, formally adopted by the IAU General Assembly in
   2000, namely MHB2000 (Mathews et al. 2002), but with the free core
   nutation (FCN - see Note 4) omitted.

4) The full MHB2000 model also contains contributions to the nutations
   in longitude and obliquity due to the free-excitation of the
   free-core-nutation during the period 1979-2000.  These FCN terms,
   which are time-dependent and unpredictable, are NOT included in the
   present function and, if required, must be independently computed.
   With the FCN corrections included, the present function delivers a
   pole which is at current epochs accurate to a few hundred
   microarcseconds.  The omission of FCN introduces further errors of
   about that size.

5) The present function provides classical nutation.  The MHB2000
   algorithm, from which it is adapted, deals also with (i) the
   offsets between the GCRS and mean poles and (ii) the adjustments in
   longitude and obliquity due to the changed precession rates.  These
   additional functions, namely frame bias and precession adjustments,
   are supported by the ERFA functions eraBi00 and eraPr00.

6) The MHB2000 algorithm also provides "total" nutations, comprising
   the arithmetic sum of the frame bias, precession adjustments,
   luni-solar nutation and planetary nutation.  These total nutations
   can be used in combination with an existing IAU 1976 precession
   implementation, such as eraPmat76, to deliver GCRS- to-true
   predictions of sub-mas accuracy at current dates.  However, there
   are three shortcomings in the MHB2000 model that must be taken into
   account if more accurate or definitive results are required (see
   Wallace 2002):

     (i) The MHB2000 total nutations are simply arithmetic sums,
         yet in reality the various components are successive Euler
         rotations.  This slight lack of rigor leads to cross terms
         that exceed 1 mas after a century.  The rigorous procedure
         is to form the GCRS-to-true rotation matrix by applying the
         bias, precession and nutation in that order.

    (ii) Although the precession adjustments are stated to be with
         respect to Lieske et al. (1977), the MHB2000 model does
         not specify which set of Euler angles are to be used and
         how the adjustments are to be applied.  The most literal
         and straightforward procedure is to adopt the 4-rotation
         epsilon_0, psi_A, omega_A, xi_A option, and to add DPSIPR
         to psi_A and DEPSPR to both omega_A and eps_A.

   (iii) The MHB2000 model predates the determination by Chapront
         et al. (2002) of a 14.6 mas displacement between the
         J2000.0 mean equinox and the origin of the ICRS frame.  It
         should, however, be noted that neglecting this displacement
         when calculating star coordinates does not lead to a
         14.6 mas change in right ascension, only a small second-
         order distortion in the pattern of the precession-nutation
         effect.

   For these reasons, the ERFA functions do not generate the "total
   nutations" directly, though they can of course easily be generated
   by calling eraBi00, eraPr00 and the present function and adding the
   results.

7) The MHB2000 model contains 41 instances where the same frequency
   appears multiple times, of which 38 are duplicates and three are
   triplicates.  To keep the present code close to the original MHB
   algorithm, this small inefficiency has not been corrected.

# References

Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
Astron.Astrophys. 387, 700

Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
Astron.Astrophys. 58, 1-16

Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.  107,
B4.  The MHB_2000 code itself was obtained on 9th September 2002 from
ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.

Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
Astron.Astrophys.Supp.Ser. 135, 111

Wallace, P.T., "Software for Implementing the IAU 2000 Resolutions",
in IERS Workshop 5.1 (2002)
"""
function nut00a(day1::AbstractFloat, day2::AbstractFloat)
    #   Interval between fundamental data J2000.0 and given date (JC.)
    Δt = ((day1 - JD2000) + day2)/(100*DAYPERYEAR)

    FACTOR_MICROARCSEC = 1 / 3.6e10

    ####    Luni-Solar Nutation
    #
    #  Fundamental (Delaunay) arguments
    #
    #  Mean anomaly of the Moon (IERS 2003).
    l = Polynomial(l0_2003A...)(Δt)
    #  Mean anomaly of the Sun (MHB 2000).
    lp = Polynomial(l1_2000A...)(Δt)
    #  Mean longitude of the Moon minus that of the ascending node (IERS 2003).
    f = Polynomial(F_2003A...)(Δt)
    #  Mean elongation of the Moon from the Sun (MHB2000).
    d = Polynomial(D_2000A...)(Δt)
    #  Mean longitude of the ascending node of the Moon (IERS 2003).
    ω = Polynomial(Ω_2003A...)(Δt)

    #  Summation of luni-solar nutation series.
    ln = ln_2000A_nutation
    la = la_2000A_nutation

    #  Convert from 0.1 μas to radians
    @inbounds begin
        arg1 = deg2rad(rem(l, ARCSECPER2PI) / 3600)
        arg2 = deg2rad(rem(lp, ARCSECPER2PI) / 3600)
        arg3 = deg2rad(rem(f, ARCSECPER2PI) / 3600)
        arg4 = deg2rad(rem(d, ARCSECPER2PI) / 3600)
        arg5 = deg2rad(rem(ω, ARCSECPER2PI) / 3600)
    end
    ϕl = Vector{Float64}(undef, size(ln, 1))
    @inbounds for i in axes(ln, 1)
        angle = ln[i,1] * arg1 + ln[i,2] * arg2 + ln[i,3] * arg3 + ln[i,4] * arg4 + ln[i,5] * arg5
        ϕl[i] = mod2pi(angle)
    end

    sum1 = sum2 = zero(eltype(la))
    @inbounds for i in axes(la, 1)
        s, c = sincos(ϕl[i])
        sum1 += (la[i,1] + la[i,2] * Δt) * s + la[i,3] * c
        sum2 += (la[i,4] + la[i,5] * Δt) * c + la[i,6] * s
    end
    δψl = deg2rad(sum1 * FACTOR_MICROARCSEC)
    δϵl = deg2rad(sum2 * FACTOR_MICROARCSEC)

    ####    Planetary Nutation
    #
    #  n.b.  The MHB2000 code computes the luni-solar and planetary
    #  nutation in different functions, using slightly different
    #  Delaunay arguments in the two cases.  This behaviour is
    #  faithfully reproduced here.  Use of the IERS 2003 expressions
    #  for both cases leads to negligible changes, well below 0.1
    #  microarcsecond.

    #  Mean anomaly of the Moon (MHB 2000).
    l00 = mod2pi(Polynomial(l0_2000A_planet...)(Δt))
    #  Mean longitude of the Moon minus that of the ascending node (MHB 2000).
    f00 = mod2pi(Polynomial(F_2000A_planet...)(Δt))
    #  Mean elongation of the Moon from the Sun (MBH 2000).
    d00 = mod2pi(Polynomial(D_2000A_planet...)(Δt))
    #  Mean longitude of the ascending node of the Moon (MHB 2000).
    ω00 = mod2pi(Polynomial(Ω_2000A_planet...)(Δt))
    #  Planetary longitudes, Mercury through Uranus (IERS 2003).
    fme = mod2pi(Polynomial(lme_2003...)(Δt))
    fve = mod2pi(Polynomial(lve_2003...)(Δt))
    fea = mod2pi(Polynomial(lea_2003...)(Δt))
    fma = mod2pi(Polynomial(lma_2003...)(Δt))
    fju = mod2pi(Polynomial(lju_2003...)(Δt))
    fsa = mod2pi(Polynomial(lsa_2003...)(Δt))
    fur = mod2pi(Polynomial(lur_2003...)(Δt))
    #  Neptune longitude (MHB 2000).
    fne = mod2pi(Polynomial(lne_2003mhb...)(Δt))
    #  General accumulated precession in longitude (IERS 2003).
    fpa = Polynomial(lge_2003...)(Δt)
    
    pn = pn_2000A_nutation
    pa = pa_2000A_nutation

    planet_args = (l00, f00, d00, ω00, fme, fve, fea, fma, fju, fsa, fur, fne, fpa)
    ϕp = Vector{Float64}(undef, size(pn, 1))
    @inbounds for i in axes(pn, 1)
        angle = zero(eltype(planet_args))
        for j in 1:length(planet_args)
            angle += pn[i,j] * planet_args[j]
        end
        ϕp[i] = mod2pi(angle)
    end

    #  Convert from 0.1 μas to radians
    sum1 = sum2 = zero(eltype(pa))
    @inbounds for i in axes(pa, 1)
        s, c = sincos(ϕp[i])
        sum1 += pa[i,1] * s + pa[i,2] * c
        sum2 += pa[i,3] * s + pa[i,4] * c
    end
    δψp = deg2rad(sum1 * FACTOR_MICROARCSEC)
    δϵp = deg2rad(sum2 * FACTOR_MICROARCSEC)

    (ψ = δψl + δψp, ϵ = δϵl + δϵp)
end

"""
    nut00b(day1::AbstractFloat, day2::AbstractFloat)

Nutation, IAU 2000B model.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `dpsi, deps` -- nutation, luni-solar + planetary (Note 2)

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
   good compromises between resolution and convenience.

2) The nutation components in longitude and obliquity are in radians
   and with respect to the equinox and ecliptic of date.  The
   obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
   value of 84381.448 arcsec.  (The errors that result from using this
   function with the IAU 2006 value of 84381.406 arcsec can be
   neglected.)

   The nutation model consists only of luni-solar terms, but includes
   also a fixed offset which compensates for certain long- period
   planetary terms (Note 7).

3) This function is an implementation of the IAU 2000B abridged
   nutation model formally adopted by the IAU General Assembly in
   2000.  The function computes the MHB_2000_SHORT luni-solar nutation
   series (Luzum 2001), but without the associated corrections for the
   precession rate adjustments and the offset between the GCRS and
   J2000.0 mean poles.

4) The full IAU 2000A (MHB2000) nutation model contains nearly 1400
   terms.  The IAU 2000B model (McCarthy & Luzum 2003) contains only
   77 terms, plus additional simplifications, yet still delivers
   results of 1 mas accuracy at present epochs.  This combination of
   accuracy and size makes the IAU 2000B abridged nutation model
   suitable for most practical applications.

   The function delivers a pole accurate to 1 mas from 1900 to 2100
   (usually better than 1 mas, very occasionally just outside 1 mas).
   The full IAU 2000A model, which is implemented in the function
   eraNut00a (q.v.), delivers considerably greater accuracy at current
   dates; however, to realize this improved accuracy, corrections for
   the essentially unpredictable free-core-nutation (FCN) must also be
   included.

5) The present function provides classical nutation.  The
   MHB_2000_SHORT algorithm, from which it is adapted, deals also with
   (i) the offsets between the GCRS and mean poles and (ii) the
   adjustments in longitude and obliquity due to the changed
   precession rates.  These additional functions, namely frame bias
   and precession adjustments, are supported by the ERFA functions
   eraBi00 and eraPr00.

6) The MHB_2000_SHORT algorithm also provides "total" nutations,
   comprising the arithmetic sum of the frame bias, precession
   adjustments, and nutation (luni-solar + planetary).  These total
   nutations can be used in combination with an existing IAU 1976
   precession implementation, such as eraPmat76, to deliver GCRS-
   to-true predictions of mas accuracy at current epochs.  However,
   for symmetry with the eraNut00a function (q.v. for the reasons),
   the ERFA functions do not generate the "total nutations" directly.
   Should they be required, they could of course easily be generated
   by calling eraBi00, eraPr00 and the present function and adding the
   results.

7) The IAU 2000B model includes "planetary bias" terms that are fixed
   in size but compensate for long-period nutations.  The amplitudes
   quoted in McCarthy & Luzum (2003), namely Dpsi = -1.5835 mas and
   Depsilon = +1.6339 mas, are optimized for the "total nutations"
   method described in Note 6.  The Luzum (2001) values used in this
   ERFA implementation, namely -0.135 mas and +0.388 mas, are
   optimized for the "rigorous" method, where frame bias, precession
   and nutation are applied separately and in that order.  During the
   interval 1995-2050, the ERFA implementation delivers a maximum
   error of 1.001 mas (not including FCN).

# References

Lieske, J.H., Lederle, T., Fricke, W., Morando, B., "Expressions for
the precession quantities based upon the IAU /1976/ system of
astronomical constants", Astron.Astrophys. 58, 1-2, 1-16. (1977)

Luzum, B., private communication, 2001 (Fortran code MHB_2000_SHORT)

McCarthy, D.D. & Luzum, B.J., "An abridged model of the
precession-nutation of the celestial pole", Cel.Mech.Dyn.Astron.  85,
37-49 (2003)

Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
Francou, G., Laskar, J., Astron.Astrophys. 282, 663-683 (1994)
"""
function nut00b(day1::AbstractFloat, day2::AbstractFloat)
    #  Interval between fundamental date J2000.0 and given date (JC).
    Δt = ((day1 - JD2000) + day2)/(100*DAYPERYEAR)

    ####    Luni-Solar Nutation
    #
    #  Fundamental (Delaunay) arguments from Simon et al. (1994).
    #
    #  Mean anomaly of the Moon.
    l = Polynomial(l0_2000B...)(Δt)
    #  Mean anomaly of the Sun.
    lp = Polynomial(l1_2000B...)(Δt)
    #  Mean longitude of the Moon minus that of the ascending node.
    f = Polynomial(F_2000B...)(Δt)
    #  Mean elongation of the Moon from the Sun.
    d = Polynomial(D_2000B...)(Δt)
    #  Mean longitude of the ascending node of the Moon.
    ω = Polynomial(Ω_2000B...)(Δt)

    #  Summation of luni-solar nutation series.
    ln = ln_2000B_nutation
    la = la_2000B_nutation

    @inbounds begin
        arg1 = deg2rad(rem(l, ARCSECPER2PI) / 3600)
        arg2 = deg2rad(rem(lp, ARCSECPER2PI) / 3600)
        arg3 = deg2rad(rem(f, ARCSECPER2PI) / 3600)
        arg4 = deg2rad(rem(d, ARCSECPER2PI) / 3600)
        arg5 = deg2rad(rem(ω, ARCSECPER2PI) / 3600)
    end

    sum1 = sum2 = zero(eltype(la))
    @inbounds for i in axes(ln, 1)
        angle = ln[i,1] * arg1 + ln[i,2] * arg2 + ln[i,3] * arg3 + ln[i,4] * arg4 + ln[i,5] * arg5
        angle = mod2pi(angle)
        s, c = sincos(angle)
        sum1 += (la[i,1] + la[i,2] * Δt) * s + la[i,3] * c
        sum2 += (la[i,4] + la[i,5] * Δt) * c + la[i,6] * s
    end

    #  Convert from 0.1 μas to radians
    DEG2RAD_FACTOR = deg2rad(1 / 3.6e10)
    δψl = sum1 * DEG2RAD_FACTOR
    δϵl = sum2 * DEG2RAD_FACTOR

    ####    In lieu of Planetary Nutation
    #
    #  Fixed offset to correct for missing terms in truncated series
    DEG2RAD_PLANET_FACTOR = deg2rad(1 / 3.6e6)
    δψp = ψ_2000B_planet * DEG2RAD_PLANET_FACTOR
    δϵp = ϵ_2000B_planet * DEG2RAD_PLANET_FACTOR

    (ψ = δψl + δψp, ϵ = δϵl + δϵp)
end

"""
    nut06a(day1::AbstractFloat, day2::AbstractFloat)

IAU 2000A nutation with adjustments to match the IAU 2006 precession.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `dpsi, deps` -- nutation, luni-solar + planetary (Note 2)

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
   good compromises between resolution and convenience.

2) The nutation components in longitude and obliquity are in radians
   and with respect to the mean equinox and ecliptic of date, IAU 2006
   precession model (Hilton et al. 2006, Capitaine et al.  2005).

3) The function first computes the IAU 2000A nutation, then applies
   adjustments for (i) the consequences of the change in obliquity
   from the IAU 1980 ecliptic to the IAU 2006 ecliptic and (ii) the
   secular variation in the Earth's dynamical form factor J2.

4) The present function provides classical nutation, complementing the
   IAU 2000 frame bias and IAU 2006 precession.  It delivers a pole
   which is at current epochs accurate to a few tens of
   microarcseconds, apart from the free core nutation.

# References

Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
Astron.Astrophys. 387, 700

Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
Astron.Astrophys. 58, 1-16

Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.  107,
B4.  The MHB_2000 code itself was obtained on 9th September 2002 from
ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.

Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
Astron.Astrophys.Supp.Ser. 135, 111

Wallace, P.T., "Software for Implementing the IAU 2000 Resolutions",
in IERS Workshop 5.1 (2002)
"""
function nut06a(day1::AbstractFloat, day2::AbstractFloat)
    #  Interval between fundamental date J2000.0 and given date (JC).
    Δt = ((day1 - JD2000) + day2)/(100*DAYPERYEAR)
    #  Obtain IAU 2000A nutation
    #  Factor correcting for secular variation of J2.
    #  Apply P03 adjustments (Wallace & Capitaine, 2006, Eqs. 5)
    ψ, ϵ = values(nut00a(day1, day2)) .* ((1.0 + j2_corr_2000*Δt) .+ (p03_2000, 0.0))
    (ψ = ψ, ϵ = ϵ)
end

"""
    nut80(day1::AbstractFloat, day2::AbstractFloat)

Nutation, IAU 1980 model.

# Input

 - `day1`   -- TT as Julian Date (Note 1)
 - `day2`   -- ... Julian Date

# Output

 - `dpsi`   -- nutation in longitude (radians)
 - `deps`   -- nutation in obliquity (radians)

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
   good compromises between resolution and convenience.

2) The nutation components are with respect to the ecliptic of date.

# References

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992), Section 3.222
(p111).
"""
function nut80(day1::AbstractFloat, day2::AbstractFloat)
    #  Interval between fundamental date J2000.0 and given date (JC).
    Δt = ((day1 - JD2000) + day2)/(100*DAYPERYEAR)

    ####    Fundamental arguments
    #
    #  
    #  The nean longitude of the Moon minus the mean longitude of
    #  the Moon's perigee.
    l = deg2rad(Polynomial(l0_1980...)(Δt)/3600.0) + 2π*rem(l0_1980t*Δt, 1.0)
    #  The mean longitude of the Sun minus the mean longitude of
    #  the Sun's perigee.
    lp = deg2rad(Polynomial(l1_1980...)(Δt)/3600.0) + 2π*rem(l1_1980t*Δt, 1.0)
    #  The mean longitude of the Moon minus the mean longitude of
    #  the Moon's node.
    f = deg2rad(Polynomial(F_1980...)(Δt)/3600.0) + 2π*rem(F_1980t*Δt, 1.0)
    #  The mean elongation of the Moon from the Sun.
    d = deg2rad(Polynomial(D_1980...)(Δt)/3600.0) + 2π*rem(D_1980t*Δt, 1.0)
    #  The mean longitude of the ascending node of the lunar orbit on
    #  the ecliptic, measured from the mean equinox of data.
    ω = deg2rad(Polynomial(Ω_1980...)(Δt)/3600.0) + 2π*rem(Ω_1980t*Δt, 1.0)

    #  Summation of luni-solar nutation series.
    ln = n_1980_nutation
    la = a_1980_nutation

    @inbounds begin
        arg1 = rem2pi(l, RoundNearest)
        arg2 = rem2pi(lp, RoundNearest)
        arg3 = rem2pi(f, RoundNearest)
        arg4 = rem2pi(d, RoundNearest)
        arg5 = rem2pi(ω, RoundNearest)
    end

    sum1 = sum2 = zero(eltype(la))
    @inbounds for i in axes(ln, 1)
        angle = ln[i,1] * arg1 + ln[i,2] * arg2 + ln[i,3] * arg3 + ln[i,4] * arg4 + ln[i,5] * arg5
        s, c = sincos(angle)
        sum1 += (la[i,1] + la[i,2] * Δt) * s
        sum2 += (la[i,3] + la[i,4] * Δt) * c
    end

    # Convert from 0.1 μas to radians
    DEG2RAD_FACTOR = deg2rad(1 / 3.6e7)
    δψl = sum1 * DEG2RAD_FACTOR
    δϵl = sum2 * DEG2RAD_FACTOR

    (ψ = δψl, ϵ = δϵl)
end

"""
    nutm80(day1::AbstractFloat, day2::AbstractFloat)

Form the matrix of nutation for a given date, IAU 1980 model.

# Input

 - `day1`  -- TDB date (Note 1)
 - `day2`  -- ... TDB date

# Output
   rmatn          double[3][3]    nutation matrix

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
   good compromises between resolution and convenience.

2) The matrix operates in the sense V(true) = rmatn * V(mean), where
   the p-vector V(true) is with respect to the true equatorial triad
   of date and the p-vector V(mean) is with respect to the mean
   equatorial triad of date.
"""
function nutm80(day1::AbstractFloat, day2::AbstractFloat)
    #  Nutation components and mean obliquity, and rotation matrix
    numat(obl80(day1, day2), values(nut80(day1, day2))...)
end


"""
    obl06(day1::AbstractFloat, day2::AbstractFloat)

Mean obliquity of the ecliptic, IAU 2006 precession model.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `ϵ`     -- obliquity of the ecliptic (radians, Note 2)

# Note

1) The TT date day1+day2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TT)=2450123.7 could be expressed in any of these ways, among
   others:

          day1          day2

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

2) The result is the angle between the ecliptic and mean equator of
   date day1+day2.

# References

Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
"""
function obl06(day1::AbstractFloat, day2::AbstractFloat)
    deg2rad(Polynomial(ϵB_2006...)(((day1-JD2000) + day2)/(100*DAYPERYEAR))/3600.0)
end

"""
    obl80(day1::AbstractFloat, day2::AbstractFloat)

Mean obliquity of the ecliptic, IAU 1980 model.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `obl`   -- obliquity of the ecliptic (radians, Note 2)

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
   good compromises between resolution and convenience.

2) The result is the angle between the ecliptic and mean equator of
   date date1+date2.

# References

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992), Expression 3.222-1
(p114).

"""
function obl80(day1::AbstractFloat, day2::AbstractFloat)
    deg2rad(Polynomial(ϵ_1980...)(((day1-JD2000) + day2)/(100*DAYPERYEAR))/3600.0)
end

"""
    p06e(day1::AbstractFloat, day2::AbstractFloat)

Precession angles, IAU 2006, equinox based.

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
   good compromises between resolution and convenience.

2) This function returns the set of equinox based angles for the
   Capitaine et al. "P03" precession theory, adopted by the IAU in
   2006.  The angles are set out in Table 1 of Hilton et al. (2006):

   eps0   epsilon_0   obliquity at J2000.0
   psia   psi_A       luni-solar precession
   oma    omega_A     inclination of equator wrt J2000.0 ecliptic
   bpa    P_A         ecliptic pole x, J2000.0 ecliptic triad
   bqa    Q_A         ecliptic pole -y, J2000.0 ecliptic triad
   pia    pi_A        angle between moving and J2000.0 ecliptics
   bpia   Pi_A        longitude of ascending node of the ecliptic
   epsa   epsilon_A   obliquity of the ecliptic
   chia   chi_A       planetary precession
   za     z_A         equatorial precession: -3rd 323 Euler angle
   zetaa  zeta_A      equatorial precession: -1st 323 Euler angle
   thetaa theta_A     equatorial precession: 2nd 323 Euler angle
   pa     p_A         general precession (n.b. see below)
   gam    gamma_J2000 J2000.0 RA difference of ecliptic poles
   phi    phi_J2000   J2000.0 codeclination of ecliptic pole
   psi    psi_J2000   longitude difference of equator poles, J2000.0

   The returned values are all radians.

   Note that the t^5 coefficient in the series for p_A from Capitaine
   et al. (2003) is incorrectly signed in Hilton et al.  (2006).

3) Hilton et al. (2006) Table 1 also contains angles that depend on
   models distinct from the P03 precession theory itself, namely the
   IAU 2000A frame bias and nutation.  The quoted polynomials are used
   in other ERFA functions:

   . eraXy06  contains the polynomial parts of the X and Y series.

   . eraS06  contains the polynomial part of the s+XY/2 series.

   . eraPfw06  implements the series for the Fukushima-Williams
     angles that are with respect to the GCRS pole (i.e. the variants
     that include frame bias).

4) The IAU resolution stipulated that the choice of parameterization
   was left to the user, and so an IAU compliant precession
   implementation can be constructed using various combinations of the
   angles returned by the present function.

5) The parameterization used by ERFA is the version of the Fukushima-
   Williams angles that refers directly to the GCRS pole.  These
   angles may be calculated by calling the function eraPfw06.  ERFA
   also supports the direct computation of the CIP GCRS X,Y by series,
   available by calling eraXy06.

6) The agreement between the different parameterizations is at the 1
   microarcsecond level in the present era.

7) When constructing a precession formulation that refers to the GCRS
   pole rather than the dynamical pole, it may (depending on the
   choice of angles) be necessary to introduce the frame bias
   explicitly.

8) It is permissible to re-use the same variable in the returned
   arguments.  The quantities are stored in the stated order.

# References

Capitaine, N., Wallace, P.T. & Chapront, J., 2003, Astron.Astrophys.,
412, 567

Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
"""
function p06e(day1::AbstractFloat, day2::AbstractFloat)
    #  Interval between fundamental date J2000.0 and given date (Julian centuries).
    Δt = ((day1 - JD2000) + day2)/(100*DAYPERYEAR)

    ####    Luni-solar precession
    #
    ψA = Polynomial(ψA_2006...)(Δt)
    #  Inclination of mean equator with respect to the J2000.0 ecliptic
    ωA = Polynomial(ωA_2006...)(Δt)
    #  Ecliptic pole x and y, J2000.0 ecliptic triad
    PA, QA = Polynomial(PA_2006...)(Δt), Polynomial(QA_2006...)(Δt)
    #  Angle between moving and J2000.0 ecliptics.
    πA = Polynomial(πA_2006...)(Δt)
    #  Longitude of the ascending node of the moving ecliptic.
    ΠA = Polynomial(ΠA_2006...)(Δt)
    #  Mean obliquity of the ecliptic.
    ϵA = Polynomial(ϵA_2006...)(Δt)

    ####    Planetary precession
    #
    χA = Polynomial(χA_2006...)(Δt)
    #  Equatorial precession: minus the first of the 323 Euler angles.
    ζA = Polynomial(ζA_2006...)(Δt)
    #  Equatorial precession: minus the second of the 323 Euler angles.
    θA = Polynomial(θA_2006...)(Δt)
    #  Equatorial precession: minus the third of the 323 Euler anlges.
    zA = Polynomial(zA_2006...)(Δt)
    #  General precession
    pA = Polynomial(pA_2006...)(Δt)
    #  Fukushima-Williams angles for precession
    γ  = Polynomial(γF_2006...)(Δt)
    ϕ  = Polynomial(ϕF_2006...)(Δt)
    ψ  = Polynomial(ψF_2006...)(Δt)

    NamedTuple{(:ϵ0, :ψA, :ωA, :Pa, :Qa, :πA, :ΠA, :ϵA, :χA, :ζA, :θA, :zA, :pA, :γ, :ϕ, :ψ)}(
        deg2rad.((ϵ0_2006, ψA, ωA, PA, QA, πA, ΠA, ϵA, χA, ζA, θA, zA, pA, γ, ϕ, ψ)./3600.0))
end

"""
    pb06(day1::AbstractFloat, day2::AbstractFloat)

This function forms three Euler angles which implement general
precession from epoch J2000.0, using the IAU 2006 model.  Frame bias
(the offset between ICRS and mean J2000.0) is included.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `ζ`    -- 1st rotation: radians cw around z
 - `θ`    -- 2nd rotation: radians ccw around y
 - `z`    -- 3rd rotation: radians cw around z

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
   good compromises between resolution and convenience.

2) The traditional accumulated precession angles zeta_A, z_A, theta_A
   cannot be obtained in the usual way, namely through polynomial
   expressions, because of the frame bias.  The latter means that two
   of the angles undergo rapid changes near this date.  They are
   instead the results of decomposing the precession-bias matrix
   obtained by using the Fukushima-Williams method, which does not
   suffer from the problem.  The decomposition returns values which
   can be used in the conventional formulation and which include frame
   bias.

3) The three angles are returned in the conventional order, which is
   not the same as the order of the corresponding Euler rotations.
   The precession-bias matrix is R_3(-z) x R_2(+theta) x R_3(-zeta).

4) Should zeta_A, z_A, theta_A angles be required that do not contain
   frame bias, they are available by calling the ERFA function
   eraP06e.
"""
function pb06(day1::AbstractFloat, day2::AbstractFloat)
    #  Precesion matrix via Fukushima-Williams angles
    r = pmat06(day1, day2)
    #  Solve for z, choosing the ±π alternative.
    x, y = -r[1,3] < 0.0 ? (r[1,3], -r[2,3]) : (-r[1,3], r[2,3])
    z = (x != 0.0 || y != 0.0) ? -atan(y, x) : 0.0
    #  De-rotate z out of the matrix
    r = Rz(z)*r
    ζ = r[2,2] != 0.0 || -r[2,1] != 0.0 ? -atan(-r[2,1], r[2,2]) : 0.0
    θ = r[3,3] != 0.0 || r[1,3] != 0.0 ? -atan(r[1,3], r[3,3]) : 0.0
    (ζ = ζ, θ = θ, z = z)
end

"""
    pfw06(day1::AbstractFloat, day2::AbstractFloat)

Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `gamb`  -- F-W angle gamma_bar (radians)
 - `phib`  -- F-W angle phi_bar (radians)
 - `psib`  -- F-W angle psi_bar (radians)
 - `epsa`  -- F-W angle epsilon_A (radians)

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
   good compromises between resolution and convenience.

2) Naming the following points:

         e = J2000.0 ecliptic pole,
         p = GCRS pole,
         E = mean ecliptic pole of date,
   and   P = mean pole of date,

   the four Fukushima-Williams angles are as follows:

      gamb = gamma_bar = epE
      phib = phi_bar = pE
      psib = psi_bar = pEP
      epsa = epsilon_A = EP

3) The matrix representing the combined effects of frame bias and
   precession is:

      PxB = R_1(-epsa).R_3(-psib).R_1(phib).R_3(gamb)

4) The matrix representing the combined effects of frame bias,
   precession and nutation is simply:

      NxPxB = R_1(-epsa-dE).R_3(-psib-dP).R_1(phib).R_3(gamb)

   where dP and dE are the nutation components with respect to the
   ecliptic of date.

# References

Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
"""
function pfw06(day1::AbstractFloat, day2::AbstractFloat)
    Δt = ((day1 - JD2000) + day2)/(100*DAYPERYEAR)

    NamedTuple{(:γ, :ϕ, :ψ, :ϵ)}(deg2rad.(
        (Polynomial(γB_2006...)(Δt), Polynomial(ϕB_2006...)(Δt),
         Polynomial(ψB_2006...)(Δt), Polynomial(ϵB_2006...)(Δt))./3600))
end

"""
    pmat00(day1::AbstractFloat, day2::AbstractFloat)

Precession matrix (including frame bias) from GCRS to a specified
date, IAU 2000 model.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `rbp`   -- bias-precession matrix (Note 2)

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
   good compromises between resolution and convenience.

2) The matrix operates in the sense V(date) = rbp * V(GCRS), where the
   p-vector V(GCRS) is with respect to the Geocentric Celestial
   Reference System (IAU, 2000) and the p-vector V(date) is with
   respect to the mean equatorial triad of the given date.

# References

IAU: Trans. International Astronomical Union, Vol. XXIVB; Proc.  24th
General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.  (2000)
"""
function pmat00(day1::AbstractFloat, day2::AbstractFloat)
    #  Obtain the required matrix (discarding others).
    bp00(day1, day2)[:rbp]
end

"""
    pmat06(day1::AbstractFloat, day2::AbstractFloat)

Precession matrix (including frame bias) from GCRS to a specified
date, IAU 2006 model.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `r`     -- bias-precession matrix (Note 2)

# Note

1) The TT date day1+day2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TT)=2450123.7 could be expressed in any of these ways, among
   others:

          day1          day2

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

2) The matrix operates in the sense V(date) = rbp * V(GCRS), where the
   p-vector V(GCRS) is with respect to the Geocentric Celestial
   Reference System (IAU, 2000) and the p-vector V(date) is with
   respect to the mean equatorial triad of the given date.

# References

Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

IAU: Trans. International Astronomical Union, Vol. XXIVB; Proc.  24th
General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.  (2000)

Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
"""
function pmat06(day1::AbstractFloat, day2::AbstractFloat)
    Δt = ((day1-JD2000) + day2)/(100*DAYPERYEAR)
    fw2m(deg2rad(1/3600).*(Polynomial(γB_2006...)(Δt),
                           Polynomial(ϕB_2006...)(Δt),
                           Polynomial(ψB_2006...)(Δt),
                           Polynomial(ϵB_2006...)(Δt))...)
end

"""
    pmat76(day1::AbstractFloat, day2::AbstractFloat)

Precession matrix from J2000.0 to a specified date, IAU 1976 model.

# Input

 - `day1`  -- ending date, TT (Note 1)
 - `day2`  -- ending date, TT (Note 1)

# Output

 - `rmatp` -- precession matrix, J2000.0 -> date1+date2

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
   good compromises between resolution and convenience.

2) The matrix operates in the sense V(date) = RMATP * V(J2000), where
   the p-vector V(J2000) is with respect to the mean equatorial triad
   of epoch J2000.0 and the p-vector V(date) is with respect to the
   mean equatorial triad of the given date.

3) Though the matrix method itself is rigorous, the precession angles
   are expressed through canonical polynomials which are valid only
   for a limited time span.  In addition, the IAU 1976 precession rate
   is known to be imperfect.  The absolute accuracy of the present
   formulation is better than 0.1 arcsec from 1960AD to 2040AD, better
   than 1 arcsec from 1640AD to 2360AD, and remains below 3 arcsec for
   the whole of the period 500BC to 3000AD.  The errors exceed 10
   arcsec outside the range 1200BC to 3900AD, exceed 100 arcsec
   outside 4200BC to 5600AD and exceed 1000 arcsec outside 6800BC to
   8200AD.

# References

Lieske, J.H., 1979, Astron.Astrophys. 73, 282.  equations (6) & (7),
p283.

Kaplan,G.H., 1981. USNO circular no. 163, pA2.
"""
function pmat76(day1::AbstractFloat, day2::AbstractFloat)
    #  Precession Euler angles, J2000.0 to specified date.
    ζ, z, θ = prec76(JD2000, 0.0, day1, day2)
    Rz(-z)Ry(θ)Rz(-ζ)
end

"""
    pn00(day1::AbstractFloat, day2::AbstractFloat, ψ::AbstractFloat, ϵ::AbstractFloat)

Precession-nutation, IAU 2000 model: a multi-purpose function,
supporting classical (equinox-based) use directly and CIO-based use
indirectly.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date
 - `dpsi,deps` -- nutation (Note 2)

# Output

 - `epsa`  -- mean obliquity (Note 3)
 - `rb`    -- frame bias matrix (Note 4)
 - `rp`    -- precession matrix (Note 5)
 - `rbp`   -- bias-precession matrix (Note 6)
 - `rn`    -- nutation matrix (Note 7)
 - `rbpn`  -- GCRS-to-true matrix (Note 8)

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
   good compromises between resolution and convenience.

2) The caller is responsible for providing the nutation components;
   they are in longitude and obliquity, in radians and are with
   respect to the equinox and ecliptic of date.  For high-accuracy
   applications, free core nutation should be included as well as any
   other relevant corrections to the position of the CIP.

3) The returned mean obliquity is consistent with the IAU 2000
   precession-nutation models.

4) The matrix rb transforms vectors from GCRS to J2000.0 mean equator
   and equinox by applying frame bias.

5) The matrix rp transforms vectors from J2000.0 mean equator and
   equinox to mean equator and equinox of date by applying precession.

6) The matrix rbp transforms vectors from GCRS to mean equator and
   equinox of date by applying frame bias then precession.  It is the
   product rp x rb.

7) The matrix rn transforms vectors from mean equator and equinox of
   date to true equator and equinox of date by applying the nutation
   (luni-solar + planetary).

8) The matrix rbpn transforms vectors from GCRS to true equator and
   equinox of date.  It is the product rn x rbp, applying frame bias,
   precession and nutation in that order.

9) It is permissible to re-use the same array in the returned
   arguments.  The arrays are filled in the order given.

# References

Capitaine, N., Chapront, J., Lambert, S. and Wallace, P., "Expressions
for the Celestial Intermediate Pole and Celestial Ephemeris Origin
consistent with the IAU 2000A precession- nutation model",
Astron.Astrophys. 400, 1145-1154 (2003)

n.b. The celestial ephemeris origin (CEO) was renamed "celestial
     intermediate origin" (CIO) by IAU 2006 Resolution 2.
"""
function pn00(day1::AbstractFloat, day2::AbstractFloat, ψ::AbstractFloat, ϵ::AbstractFloat)
    #  IAU 2000 precession-rate adjustments
    ψpr, ϵpr = pr00(day1, day2)
    #  Mean obliquity, consistent with IAU 2000 precession-nutation
    ϵA = obl80(day1, day2) + ϵpr
    #  Frame bias and precession matrices and their product.
    rb, rp, rbp = bp00(day1, day2)
    #  Nutation matrix
    rn = numat(ϵA, ψ, ϵ)
    #  Bias-precession-nutation matrix (classical)
    (ϵA = ϵA, rb = rb, rp = rp, rbp = rbp, rn = rn, rbpn = rn*rbp)
end

"""
    pn00a(day1::AbstractFloat, day2::AbstractFloat)

Precession-nutation, IAU 2000A model: a multi-purpose function,
supporting classical (equinox-based) use directly and CIO-based use
indirectly.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2   -- ... as Julian Date

# Output

 - `dpsi, deps` -- nutation (Note 2)
 - `epsa`  -- mean obliquity (Note 3)
 - `rb`    -- frame bias matrix (Note 4)
 - `rp`    -- precession matrix (Note 5)
 - `rbp`   -- bias-precession matrix (Note 6)
 - `rn`    -- nutation matrix (Note 7)
 - `rbpn`  -- GCRS-to-true matrix (Notes 8,9)

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
   good compromises between resolution and convenience.

2) The nutation components (luni-solar + planetary, IAU 2000A) in
   longitude and obliquity are in radians and with respect to the
   equinox and ecliptic of date.  Free core nutation is omitted; for
   the utmost accuracy, use the eraPn00 function, where the nutation
   components are caller-specified.  For faster but slightly less
   accurate results, use the eraPn00b function.

3) The mean obliquity is consistent with the IAU 2000 precession.

4) The matrix rb transforms vectors from GCRS to J2000.0 mean equator
   and equinox by applying frame bias.

5) The matrix rp transforms vectors from J2000.0 mean equator and
   equinox to mean equator and equinox of date by applying precession.

6) The matrix rbp transforms vectors from GCRS to mean equator and
   equinox of date by applying frame bias then precession.  It is the
   product rp x rb.

7) The matrix rn transforms vectors from mean equator and equinox of
   date to true equator and equinox of date by applying the nutation
   (luni-solar + planetary).

8) The matrix rbpn transforms vectors from GCRS to true equator and
   equinox of date.  It is the product rn x rbp, applying frame bias,
   precession and nutation in that order.

9) The X,Y,Z coordinates of the IAU 2000A Celestial Intermediate Pole
   are elements (3,1-3) of the GCRS-to-true matrix, i.e. rbpn[2][0-2].

10) It is permissible to re-use the same array in the returned
    arguments.  The arrays are filled in the stated order.

# References

Capitaine, N., Chapront, J., Lambert, S. and Wallace, P., "Expressions
for the Celestial Intermediate Pole and Celestial Ephemeris Origin
consistent with the IAU 2000A precession- nutation model",
Astron.Astrophys. 400, 1145-1154 (2003)

n.b. The celestial ephemeris origin (CEO) was renamed "celestial
     intermediate origin" (CIO) by IAU 2006 Resolution 2.
"""
function pn00a(day1::AbstractFloat, day2::AbstractFloat)
    ψ, ϵ = nut00a(day1, day2)
    ϵA, rb, rp, rbp, rn, rbpn = pn00(day1, day2, ψ, ϵ)
    (ψ = ψ, ϵ = ϵ, ϵA = ϵA, rb = rb, rp = rp, rbp = rbp, rn = rn, rbpn = rbpn)
end

"""
    pn00b(day1::AbstractFloat, day2::AbstractFloat)

Precession-nutation, IAU 2000B model: a multi-purpose function,
supporting classical (equinox-based) use directly and CIO-based use
indirectly.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `dpsi, deps` -- nutation (Note 2)
 - `epsa`  -- mean obliquity (Note 3)
 - `rb`    -- frame bias matrix (Note 4)
 - `rp`    -- precession matrix (Note 5)
 - `rbp`   -- bias-precession matrix (Note 6)
 - `rn`    -- nutation matrix (Note 7)
 - `rbpn`  -- GCRS-to-true matrix (Notes 8,9)

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
   good compromises between resolution and convenience.

2) The nutation components (luni-solar + planetary, IAU 2000B) in
   longitude and obliquity are in radians and with respect to the
   equinox and ecliptic of date.  For more accurate results, but at
   the cost of increased computation, use the eraPn00a function.  For
   the utmost accuracy, use the eraPn00 function, where the nutation
   components are caller-specified.

3) The mean obliquity is consistent with the IAU 2000 precession.

4) The matrix rb transforms vectors from GCRS to J2000.0 mean equator
   and equinox by applying frame bias.

5) The matrix rp transforms vectors from J2000.0 mean equator and
   equinox to mean equator and equinox of date by applying precession.

6) The matrix rbp transforms vectors from GCRS to mean equator and
   equinox of date by applying frame bias then precession.  It is the
   product rp x rb.

7) The matrix rn transforms vectors from mean equator and equinox of
   date to true equator and equinox of date by applying the nutation
   (luni-solar + planetary).

8) The matrix rbpn transforms vectors from GCRS to true equator and
   equinox of date.  It is the product rn x rbp, applying frame bias,
   precession and nutation in that order.

9) The X,Y,Z coordinates of the IAU 2000B Celestial Intermediate Pole
   are elements (3,1-3) of the GCRS-to-true matrix, i.e. rbpn[2][0-2].

10) It is permissible to re-use the same array in the returned
    arguments.  The arrays are filled in the stated order.

# References

Capitaine, N., Chapront, J., Lambert, S. and Wallace, P., "Expressions
for the Celestial Intermediate Pole and Celestial Ephemeris Origin
consistent with the IAU 2000A precession- nutation model",
Astron.Astrophys. 400, 1145-1154 (2003).

n.b. The celestial ephemeris origin (CEO) was renamed "celestial
     intermediate origin" (CIO) by IAU 2006 Resolution 2.
"""
function pn00b(day1::AbstractFloat, day2::AbstractFloat)
    ψ, ϵ = nut00b(day1, day2)
    ϵA, rb, rp, rbp, rn, rbpn = pn00(day1, day2, ψ, ϵ)
    (ψ = ψ, ϵ = ϵ, ϵA = ϵA, rb = rb, rp = rp, rbp = rbp, rn = rn, rbpn = rbpn)
end

"""
    pn06(day1::AbstractFloat, day2::AbstractFloat, δψ::AbstractFloat, δϵ::AbstractFloat)

Precession-nutation, IAU 2006 model: a multi-purpose function,
supporting classical (equinox-based) use directly and CIO-based use
indirectly.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date
 - `dpsi`  -- nutation (Note 2)
 - `deps`  -- nutation (Note 2)

# Output

 - `epsa`  -- mean obliquity (Note 3)
 - `rb`    -- frame bias matrix (Note 4)
 - `rp`    -- precession matrix (Note 5)
 - `rbp`   -- bias-precession matrix (Note 6)
 - `rn`    -- nutation matrix (Note 7)
 - `rbpn`  -- GCRS-to-true matrix (Notes 8,9)

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
   good compromises between resolution and convenience.

2) The caller is responsible for providing the nutation components;
   they are in longitude and obliquity, in radians and are with
   respect to the equinox and ecliptic of date.  For high-accuracy
   applications, free core nutation should be included as well as any
   other relevant corrections to the position of the CIP.

3) The returned mean obliquity is consistent with the IAU 2006
   precession.

4) The matrix rb transforms vectors from GCRS to J2000.0 mean equator
   and equinox by applying frame bias.

5) The matrix rp transforms vectors from J2000.0 mean equator and
   equinox to mean equator and equinox of date by applying precession.

6) The matrix rbp transforms vectors from GCRS to mean equator and
   equinox of date by applying frame bias then precession.  It is the
   product rp x rb.

7) The matrix rn transforms vectors from mean equator and equinox of
   date to true equator and equinox of date by applying the nutation
   (luni-solar + planetary).

8) The matrix rbpn transforms vectors from GCRS to true equator and
   equinox of date.  It is the product rn x rbp, applying frame bias,
   precession and nutation in that order.

9) The X,Y,Z coordinates of the Celestial Intermediate Pole are
   elements (3,1-3) of the GCRS-to-true matrix, i.e. rbpn[2][0-2].

10) It is permissible to re-use the same array in the returned
    arguments.  The arrays are filled in the stated order.

# References

Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
"""
function pn06(day1::AbstractFloat, day2::AbstractFloat, δψ::AbstractFloat, δϵ::AbstractFloat)
    #  Bias-precession Fukushima-Williams angle of J2000.0 = frame bias
    #  and matrix
    rb = fw2m(pfw06(MJD0, MJD00)...)
    #  Bias-precession Fukushima-Williams angles of date.
    γb, ϕb, ψb, ϵb = pfw06(day1, day2)
    rbp = fw2m(γb, ϕb, ψb, ϵb)
    #  Solve for precession matrix
    rp = rbp*rb'
    #  Equinox-based bias-precession-nutation matrix
    rbpn = fw2m(γb, ϕb, ψb + δψ, ϵb + δϵ)
    #  Solve for nutation matrix
    rn = rbpn*rbp'
    (ϵb = ϵb, rb = rb, rp = rp, rbp = rbp, rn = rn, rbpn = rbpn)
end

"""
    pn06a(day1::AbstractFloat, day2::AbstractFloat)

Precession-nutation, IAU 2006/2000A models: a multi-purpose function,
supporting classical (equinox-based) use directly and CIO-based use
indirectly.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `dpsi, deps` -- nutation (Note 2)
 - `epsa`  -- mean obliquity (Note 3)
 - `rb`    -- frame bias matrix (Note 4)
 - `rp`    -- precession matrix (Note 5)
 - `rbp`   -- bias-precession matrix (Note 6)
 - `rn`    -- nutation matrix (Note 7)
 - `rbpn`  -- GCRS-to-true matrix (Notes 8,9)

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
   good compromises between resolution and convenience.

2) The nutation components (luni-solar + planetary, IAU 2000A) in
   longitude and obliquity are in radians and with respect to the
   equinox and ecliptic of date.  Free core nutation is omitted; for
   the utmost accuracy, use the eraPn06 function, where the nutation
   components are caller-specified.

3) The mean obliquity is consistent with the IAU 2006 precession.

4) The matrix rb transforms vectors from GCRS to mean J2000.0 by
   applying frame bias.

5) The matrix rp transforms vectors from mean J2000.0 to mean of date
   by applying precession.

6) The matrix rbp transforms vectors from GCRS to mean of date by
   applying frame bias then precession.  It is the product rp x rb.

7) The matrix rn transforms vectors from mean of date to true of date
   by applying the nutation (luni-solar + planetary).

8) The matrix rbpn transforms vectors from GCRS to true of date
   (CIP/equinox).  It is the product rn x rbp, applying frame bias,
   precession and nutation in that order.

9) The X,Y,Z coordinates of the IAU 2006/2000A Celestial Intermediate
   Pole are elements (3,1-3) of the GCRS-to-true matrix,
   i.e. rbpn[2][0-2].

10) It is permissible to re-use the same array in the returned
    arguments.  The arrays are filled in the stated order.

# References

Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
"""
function pn06a(day1::AbstractFloat, day2::AbstractFloat)
    ψ, ϵ = nut06a(day1, day2)
    ϵA, rb, rp, rbp, rn, rbpn = pn06(day1, day2, ψ, ϵ)
    (ψ = ψ, ϵ = ϵ, ϵA = ϵA, rb = rb, rp = rp, rbp = rbp, rn = rn, rbpn = rbpn)
end

"""
    pnm00a(day1::AbstractFloat, day2::AbstractFloat)

Form the matrix of precession-nutation for a given date (including
frame bias), equinox based, IAU 2000A model.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `rbpn`  -- bias-precession-nutation matrix (Note 2)

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
   good compromises between resolution and convenience.

2) The matrix operates in the sense V(date) = rbpn * V(GCRS), where
   the p-vector V(date) is with respect to the true equatorial triad
   of date date1+date2 and the p-vector V(GCRS) is with respect to the
   Geocentric Celestial Reference System (IAU, 2000).

3) A faster, but slightly less accurate, result (about 1 mas) can be
   obtained by using instead the eraPnm00b function.

# References

IAU: Trans. International Astronomical Union, Vol. XXIVB; Proc.  24th
General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.  (2000)

"""
function pnm00a(day1::AbstractFloat, day2::AbstractFloat)
    pn00a(day1, day2)[:rbpn]
end

"""
    pnm00b(day1::AbstractFloat, day2::AbstractFloat)

Form the matrix of precession-nutation for a given date (including
frame bias), equinox-based, IAU 2000B model.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `rbpn`  -- bias-precession-nutation matrix (Note 2)

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
   good compromises between resolution and convenience.

2) The matrix operates in the sense V(date) = rbpn * V(GCRS), where
   the p-vector V(date) is with respect to the true equatorial triad
   of date date1+date2 and the p-vector V(GCRS) is with respect to the
   Geocentric Celestial Reference System (IAU, 2000).

3) The present function is faster, but slightly less accurate (about 1
   mas), than the eraPnm00a function.

# References

IAU: Trans. International Astronomical Union, Vol. XXIVB; Proc.  24th
General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.  (2000)
"""
function pnm00b(day1::AbstractFloat, day2::AbstractFloat)
    pn00b(day1, day2)[:rbpn]
end

"""
    pnm06a(day1::AbstractFloat, day2::AbstractFloat)

Form the matrix of precession-nutation for a given date (including
frame bias), equinox based, IAU 2006 precession and IAU 2000A nutation
models.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `rbpn`  -- bias-precession-nutation matrix (Note 2)

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
   good compromises between resolution and convenience.

2) The matrix operates in the sense V(date) = rbpn * V(GCRS), where
   the p-vector V(date) is with respect to the true equatorial triad
   of date date1+date2 and the p-vector V(GCRS) is with respect to the
   Geocentric Celestial Reference System (IAU, 2000).

# References

Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855.
"""
function pnm06a(day1::AbstractFloat, day2::AbstractFloat)
    #  Fukashima-Williams angles for frame bias and precession
    γB, ϕB, ψB, ϵA = pfw06(day1, day2)
    #  Nutation components
    δϕ, δϵ = nut06a(day1, day2)
    #  Equinox based nutation x precesion x bias matrix
    fw2m(γB, ϕB, ψB + δϕ, ϵA + δϵ)
end

"""
    pnm80(day1::AbstractFloat, day2::AbstractFloat)

Form the matrix of precession/nutation for a given date, IAU 1976
precession model, IAU 1980 nutation model.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `rmatpn` -- combined precession/nutation matrix

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
   good compromises between resolution and convenience.

2) The matrix operates in the sense V(date) = rmatpn * V(J2000), where
   the p-vector V(date) is with respect to the true equatorial triad
   of date date1+date2 and the p-vector V(J2000) is with respect to
   the mean equatorial triad of epoch J2000.0.

# References

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992), Section 3.3 (p145).
"""
function pnm80(day1::AbstractFloat, day2::AbstractFloat)
    #  Precession and nutation matrices J2000.0 to date
    nutm80(day1, day2)*pmat76(day1, day2)
end

"""
    pom00(x::AbstractFloat, y::AbstractFloat, s::AbstractFloat)

Form the matrix of polar motion for a given date, IAU 2000.

# Input

 - `xp, yp` -- coordinates of the pole (radians, Note 1)
 - `sp`     -- the TIO locator s' (radians, Note 2)

# Output

 - `rpom`  -- polar-motion matrix (Note 3)

# Note

1) The arguments xp and yp are the coordinates (in radians) of the
   Celestial Intermediate Pole with respect to the International
   Terrestrial Reference System (see IERS Conventions 2003), measured
   along the meridians 0 and 90 deg west respectively.

2) The argument sp is the TIO locator s', in radians, which positions
   the Terrestrial Intermediate Origin on the equator.  It is obtained
   from polar motion observations by numerical integration, and so is
   in essence unpredictable.  However, it is dominated by a secular
   drift of about 47 microarcseconds per century, and so can be taken
   into account by using s' = -47*t, where t is centuries since
   J2000.0.  The function eraSp00 implements this approximation.

3) The matrix operates in the sense V(TRS) = rpom * V(CIP), meaning
   that it is the final rotation when computing the pointing direction
   to a celestial source.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
pom00(x::AbstractFloat, y::AbstractFloat, s::AbstractFloat) = Rx(-y)Ry(-x)Rz(s)

"""
    pr00(day1::AbstractFloat, day2::AbstractFloat)

Precession-rate part of the IAU 2000 precession-nutation models (part
of MHB2000).

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `dpsipr, depspr` -- precession corrections (Notes 2,3)

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
   good compromises between resolution and convenience.

2) The precession adjustments are expressed as "nutation components",
   corrections in longitude and obliquity with respect to the J2000.0
   equinox and ecliptic.

3) Although the precession adjustments are stated to be with respect
   to Lieske et al. (1977), the MHB2000 model does not specify which
   set of Euler angles are to be used and how the adjustments are to
   be applied.  The most literal and straightforward procedure is to
   adopt the 4-rotation epsilon_0, psi_A, omega_A, xi_A option, and to
   add dpsipr to psi_A and depspr to both omega_A and eps_A.

4) This is an implementation of one aspect of the IAU 2000A nutation
   model, formally adopted by the IAU General Assembly in 2000, namely
   MHB2000 (Mathews et al. 2002).

# References

Lieske, J.H., Lederle, T., Fricke, W. & Morando, B., "Expressions for
the precession quantities based upon the IAU (1976) System of
Astronomical Constants", Astron.Astrophys., 58, 1-16 (1977)

Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation and
precession New nutation series for nonrigid Earth and insights into
the Earth's interior", J.Geophys.Res., 107, B4, 2002.  The MHB2000
code itself was obtained on 9th September 2002 from
ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.

Wallace, P.T., "Software for Implementing the IAU 2000 Resolutions",
in IERS Workshop 5.1 (2002).
"""
function pr00(day1::AbstractFloat, day2::AbstractFloat)
    Δt = ((day1 - JD2000) + day2)/(100*DAYPERYEAR)
    #  Precession and obliquity corrections (radians/century)/
    #  Precession rate contributions with respect to IAU 1976/1980
    ψ, ϵ = deg2rad.((ψ_corr_2000, ϵ_corr_2000).*Δt./3600)
    (ψ = ψ, ϵ = ϵ)
end

"""
    prec76(day11::AbstractFloat, day12::AbstractFloat, day21::AbstractFloat, day22::AbstractFloat)

IAU 1976 precession model.

This function forms the three Euler angles which implement general
precession between two dates, using the IAU 1976 model (as for the FK5
catalog).

# Input

 - `day01` -- TDB starting date (Note 1)
 - `day02` -- ... starting date
 - `day11` -- TDB ending date (Note 1)
 - `day12` -- ... starting date

# Output

 - `zeta`  -- 1st rotation: radians cw around z
 - `z`     -- 3rd rotation: radians cw around z
 - `theta` -- 2nd rotation: radians ccw around y

# Note

1) The dates date01+date02 and date11+date12 are Julian Dates,
   apportioned in any convenient way between the arguments daten1 and
   daten2.  For example, JD(TDB)=2450123.7 could be expressed in any
   of these ways, among others:

         daten1        daten2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 method is best matched to the way the
   argument is handled internally and will deliver the optimum optimum
   resolution.  The MJD method and the date & time methods are both
   good compromises between resolution and convenience.  The two dates
   may be expressed using different methods, but at the risk of losing
   some resolution.

2) The accumulated precession angles zeta, z, theta are expressed
   through canonical polynomials which are valid only for a limited
   time span.  In addition, the IAU 1976 precession rate is known to
   be imperfect.  The absolute accuracy of the present formulation is
   better than 0.1 arcsec from 1960AD to 2040AD, better than 1 arcsec
   from 1640AD to 2360AD, and remains below 3 arcsec for the whole of
   the period 500BC to 3000AD.  The errors exceed 10 arcsec outside
   the range 1200BC to 3900AD, exceed 100 arcsec outside 4200BC to
   5600AD and exceed 1000 arcsec outside 6800BC to 8200AD.

3) The three angles are returned in the conventional order, which is
   not the same as the order of the corresponding Euler rotations.
   The precession matrix is R_3(-z) x R_2(+theta) x R_3(-zeta).

# References

Lieske, J.H., 1979, Astron.Astrophys. 73, 282, equations (6) & (7),
p283.
"""
function prec76(day11::AbstractFloat, day12::AbstractFloat, day21::AbstractFloat, day22::AbstractFloat)
    #  Interval between fundamental epoch J2000.0 and start date (Julian centuries).
    t0 = ((day11 - JD2000) + day12)/(100*DAYPERYEAR)
    #  Interval over which precession required (Julian centuries).
    Δt = ((day21 - day11) + (day22 - day12))/(100*DAYPERYEAR)
    #  Euler angles.
    wt, θt = Polynomial(ζT_1976...)(t0), Polynomial(θT_1976...)(t0)
    ζ = Polynomial(SVector(0., wt, Polynomial(ζA_1976[1:2]...)(t0), ζA_1976[3])...)(Δt)
    z = Polynomial(SVector(0., wt, Polynomial(zA_1976[1:2]...)(t0), zA_1976[3])...)(Δt)
    θ = Polynomial(SVector(0., θt, Polynomial(θA_1976[1:2]...)(t0), θA_1976[3])...)(Δt)
    ζ, z, θ = deg2rad.((ζ, z, θ)./3600.0)
    (ζ = ζ, z = z, θ = θ)
end

"""
    s00(day1::AbstractFloat, day2::AbstractFloat, x::AbstractFloat, y::AbstractFloat)

The CIO locator s, positioning the Celestial Intermediate Origin on
the equator of the Celestial Intermediate Pole, given the CIP's X,Y
coordinates.  Compatible with IAU 2000A precession-nutation.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date
 - `x, y`  -- CIP coordinates (Note 3)

# Output

 - `s`     -- the CIO locator s in radians (Note 2)

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
   good compromises between resolution and convenience.

2) The CIO locator s is the difference between the right ascensions of
   the same point in two systems: the two systems are the GCRS and the
   CIP,CIO, and the point is the ascending node of the CIP equator.
   The quantity s remains below 0.1 arcsecond throughout 1900-2100.

3) The series used to compute s is in fact for s+XY/2, where X and Y
   are the x and y components of the CIP unit vector; this series is
   more compact than a direct series for s would be.  This function
   requires X,Y to be supplied by the caller, who is responsible for
   providing values that are consistent with the supplied date.

4) The model is consistent with the IAU 2000A precession-nutation.

# References

Capitaine, N., Chapront, J., Lambert, S. and Wallace, P., "Expressions
for the Celestial Intermediate Pole and Celestial Ephemeris Origin
consistent with the IAU 2000A precession- nutation model",
Astron.Astrophys. 400, 1145-1154 (2003)

n.b. The celestial ephemeris origin (CEO) was renamed "celestial
    intermediate origin" (CIO) by IAU 2006 Resolution 2.

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function s00(day1::AbstractFloat, day2::AbstractFloat, x::AbstractFloat, y::AbstractFloat)
    Δt = ((day1 - JD2000) + day2)/(100*DAYPERYEAR)

    #  Fundamental Arguments (IERS Conventions 2003)
    @inline ϕ = SVector(
        #  Mean anomaly of the Moon
        fal03(Δt),
        #  Mean anomaly of the Sun
        falp03(Δt),
        #  Mean longitude of the Moon minus that of the ascending node
        faf03(Δt),
        #  Mean elongation of the Moon from the Sun
        fad03(Δt),
        #  Mean longitude of the ascending node of the Moon
        faom03(Δt),
        #  Mean longitude of Venus
        fave03(Δt),
        #  Mean longitude of Earth
        fae03(Δt),
        #  General precession in longitude
        fapa03(Δt))

    ϕ0 = ϕ0_2000As*ϕ
    a0 = a0_2000As
    ϕ1 = ϕ1_2000As*ϕ
    a1 = a1_2000As
    ϕ2 = ϕ2_2000As*ϕ
    a2 = a2_2000As
    ϕ3 = ϕ3_2000As*ϕ
    a3 = a3_2000As
    ϕ4 = ϕ4_2000As*ϕ
    a4 = a4_2000As

    sum0 = sum1 = sum2 = sum3 = sum4 = zero(eltype(a0))

    @inbounds for i in axes(a0, 1)
        s, c = sincos(ϕ0[i])
        sum0 += a0[i,1] * s + a0[i,2] * c
    end
    @inbounds for i in axes(a1, 1)
        s, c = sincos(ϕ1[i])
        sum1 += a1[i,1] * s + a1[i,2] * c
    end
    @inbounds for i in axes(a2, 1)
        s, c = sincos(ϕ2[i])
        sum2 += a2[i,1] * s + a2[i,2] * c
    end
    @inbounds for i in axes(a3, 1)
        s, c = sincos(ϕ3[i])
        sum3 += a3[i,1] * s + a3[i,2] * c
    end
    @inbounds for i in axes(a4, 1)
        s, c = sincos(ϕ4[i])
        sum4 += a4[i,1] * s + a4[i,2] * c
    end

    corrections = SVector(sum0, sum1, sum2, sum3, sum4, 0.0)
    deg2rad(Polynomial((sp_2000A .+ corrections)...)(Δt) / 3600) - x*y/2.0
end

"""
    s00a(day1::AbstractFloat, day2::AbstractFloat)

The CIO locator s, positioning the Celestial Intermediate Origin on
the equator of the Celestial Intermediate Pole, using the IAU 2000A
precession-nutation model.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `s`     -- the CIO locator s in radians (Note 2)

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
   good compromises between resolution and convenience.

2) The CIO locator s is the difference between the right ascensions of
   the same point in two systems.  The two systems are the GCRS and
   the CIP,CIO, and the point is the ascending node of the CIP
   equator.  The CIO locator s remains a small fraction of 1 arcsecond
   throughout 1900-2100.

3) The series used to compute s is in fact for s+XY/2, where X and Y
   are the x and y components of the CIP unit vector; this series is
   more compact than a direct series for s would be.  The present
   function uses the full IAU 2000A nutation model when predicting the
   CIP position.  Faster results, with no significant loss of
   accuracy, can be obtained via the function eraS00b, which uses
   instead the IAU 2000B truncated model.

# References

Capitaine, N., Chapront, J., Lambert, S. and Wallace, P., "Expressions
for the Celestial Intermediate Pole and Celestial Ephemeris Origin
consistent with the IAU 2000A precession- nutation model",
Astron.Astrophys. 400, 1145-1154 (2003)

n.b. The celestial ephemeris origin (CEO) was renamed "celestial
     intermediate origin" (CIO) by IAU 2006 Resolution 2.

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function s00a(day1::AbstractFloat, day2::AbstractFloat)
    #  Bias-precession-nutation matrix (IAU 2000A), extract the CIP
    #  coordinates, and compute the CIO locator s.
    s00(day1, day2, bpn2xy(pnm00a(day1, day2))...)
end

"""
    s00b(day1::AbstractFloat, day2::AbstractFloat)

The CIO locator s, positioning the Celestial Intermediate Origin on
the equator of the Celestial Intermediate Pole, using the IAU 2000B
precession-nutation model.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - 's'     -- the CIO locator s in radians (Note 2)

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
   good compromises between resolution and convenience.

2) The CIO locator s is the difference between the right ascensions of
   the same point in two systems.  The two systems are the GCRS and
   the CIP,CIO, and the point is the ascending node of the CIP
   equator.  The CIO locator s remains a small fraction of 1 arcsecond
   throughout 1900-2100.

3) The series used to compute s is in fact for s+XY/2, where X and Y
   are the x and y components of the CIP unit vector; this series is
   more compact than a direct series for s would be.  The present
   function uses the IAU 2000B truncated nutation model when
   predicting the CIP position.  The function eraS00a uses instead the
   full IAU 2000A model, but with no significant increase in accuracy
   and at some cost in speed.

# References

Capitaine, N., Chapront, J., Lambert, S. and Wallace, P., "Expressions
for the Celestial Intermediate Pole and Celestial Ephemeris Origin
consistent with the IAU 2000A precession- nutation model",
Astron.Astrophys. 400, 1145-1154 (2003)

n.b. The celestial ephemeris origin (CEO) was renamed "celestial
     intermediate origin" (CIO) by IAU 2006 Resolution 2.

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function s00b(day1::AbstractFloat, day2::AbstractFloat)
    #  Bias-precession-nutation matrix (IAU 2000A), extract the CIP
    #  coordinates, and compute the CIO locator s.
    s00(day1, day2, bpn2xy(pnm00b(day1, day2))...)
end

"""
    s06(day1::AbstractFloat, day2::AbstractFloat, x::AbstractFloat, y::AbstractFloat)

The CIO locator s, positioning the Celestial Intermediate Origin on
the equator of the Celestial Intermediate Pole, given the CIP's X,Y
coordinates.  Compatible with IAU 2006/2000A precession-nutation.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date
 - `x, y`  -- CIP coordinates (Note 3)

# Output

 - `s`     -- the CIO locator s in radians (Note 2)

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
   good compromises between resolution and convenience.

2) The CIO locator s is the difference between the right ascensions of
   the same point in two systems: the two systems are the GCRS and the
   CIP,CIO, and the point is the ascending node of the CIP equator.
   The quantity s remains below 0.1 arcsecond throughout 1900-2100.

3) The series used to compute s is in fact for s+XY/2, where X and Y
   are the x and y components of the CIP unit vector; this series is
   more compact than a direct series for s would be.  This function
   requires X,Y to be supplied by the caller, who is responsible for
   providing values that are consistent with the supplied date.

4) The model is consistent with the "P03" precession (Capitaine et
   al. 2003), adopted by IAU 2006 Resolution 1, 2006, and the IAU
   2000A nutation (with P03 adjustments).

# References

Capitaine, N., Wallace, P.T. & Chapront, J., 2003, Astron.
Astrophys. 432, 355

McCarthy, D.D., Petit, G. (eds.) 2004, IERS Conventions (2003), IERS
Technical Note No. 32, BKG
"""
function s06(day1::AbstractFloat, day2::AbstractFloat, x::AbstractFloat, y::AbstractFloat)
    Δt = ((day1 - JD2000) + day2)/(100*DAYPERYEAR)

    #  Fundamental Arguments (IERS Conventions 2003)
    @inline ϕ = SVector(
        #  Mean anomaly of the Moon
        fal03(Δt),
        #  Mean anomaly of the Sun
        falp03(Δt),
        #  Mean longitude of the Moon minus that of the ascending node
        faf03(Δt),
        #  Mean elongation of the Moon from the Sun
        fad03(Δt),
        #  Mean longitude of the ascending node of the Moon
        faom03(Δt),
        #  Mean longitude of Venus
        fave03(Δt),
        #  Mean longitude of Earth
        fae03(Δt),
        #  General precession in longitude
        fapa03(Δt))

    # a2_2006_equinox
    ϕ0 = ϕ0_2006_equinox*ϕ
    a0 = a0_2006_equinox
    ϕ1 = ϕ1_2006_equinox*ϕ
    a1 = a1_2006_equinox
    ϕ2 = ϕ2_2006_equinox*ϕ
    a2 = a2_2006_equinox
    ϕ3 = ϕ3_2006_equinox*ϕ
    a3 = a3_2006_equinox
    ϕ4 = ϕ4_2006_equinox*ϕ
    a4 = a4_2006_equinox

    sum0 = sum1 = sum2 = sum3 = sum4 = zero(eltype(a0))

    @inbounds for i in axes(a0, 1)
        s, c = sincos(ϕ0[i])
        sum0 += a0[i,1] * s + a0[i,2] * c
    end
    @inbounds for i in axes(a1, 1)
        s, c = sincos(ϕ1[i])
        sum1 += a1[i,1] * s + a1[i,2] * c
    end
    @inbounds for i in axes(a2, 1)
        s, c = sincos(ϕ2[i])
        sum2 += a2[i,1] * s + a2[i,2] * c
    end
    @inbounds for i in axes(a3, 1)
        s, c = sincos(ϕ3[i])
        sum3 += a3[i,1] * s + a3[i,2] * c
    end
    @inbounds for i in axes(a4, 1)
        s, c = sincos(ϕ4[i])
        sum4 += a4[i,1] * s + a4[i,2] * c
    end

    corrections = SVector(sum0, sum1, sum2, sum3, sum4, 0.0)
    deg2rad(Polynomial((cio_s_2006 .+ corrections)...)(Δt) / 3600) - x*y/2.0
end

"""
    s06a(day1::AbstractFloat, day2::AbstractFloat)

The CIO locator s, positioning the Celestial Intermediate Origin on
the equator of the Celestial Intermediate Pole, using the IAU 2006
precession and IAU 2000A nutation models.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `s`     -- the CIO locator s in radians (Note 2)

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
   good compromises between resolution and convenience.

2) The CIO locator s is the difference between the right ascensions of
   the same point in two systems.  The two systems are the GCRS and
   the CIP,CIO, and the point is the ascending node of the CIP
   equator.  The CIO locator s remains a small fraction of 1 arcsecond
   throughout 1900-2100.

3) The series used to compute s is in fact for s+XY/2, where X and Y
   are the x and y components of the CIP unit vector; this series is
   more compact than a direct series for s would be.  The present
   function uses the full IAU 2000A nutation model when predicting the
   CIP position.

# References

Capitaine, N., Chapront, J., Lambert, S. and Wallace, P., "Expressions
for the Celestial Intermediate Pole and Celestial Ephemeris Origin
consistent with the IAU 2000A precession- nutation model",
Astron.Astrophys. 400, 1145-1154 (2003)

n.b. The celestial ephemeris origin (CEO) was renamed "celestial
     intermediate origin" (CIO) by IAU 2006 Resolution 2.

Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003), IERS
Technical Note No. 32, BKG

Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
"""
function s06a(day1::AbstractFloat, day2::AbstractFloat)
    #  Bias-precession-nutation matrix (IAU 2006/2000A), extract the CIP
    #  coordinates, and compute the CIO locator s.
    s06(day1, day2, bpn2xy(pnm06a(day1, day2))...)
end

"""
    sp00(day1::AbstractFloat, day2::AbstractFloat)

The TIO locator s', positioning the Terrestrial Intermediate Origin on
the equator of the Celestial Intermediate Pole.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `s`     -- the TIO locator s' in radians (Note 2)

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
   good compromises between resolution and convenience.

2) The TIO locator s' is obtained from polar motion observations by
   numerical integration, and so is in essence unpredictable.
   However, it is dominated by a secular drift of about 47
   microarcseconds per century, which is the approximation evaluated
   by the present function.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function sp00(day1::AbstractFloat, day2::AbstractFloat)
    deg2rad(tio_2000*((day1 - JD2000) + day2)/(100*DAYPERYEAR)/3600.0)
end

"""
    xy06(day1::AbstractFloat, day2::AbstractFloat)

X,Y coordinates of celestial intermediate pole from series based on
IAU 2006 precession and IAU 2000A nutation.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `x, y`  -- CIP X,Y coordinates (Note 2)

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
   good compromises between resolution and convenience.

2) The X,Y coordinates are those of the unit vector towards the
   celestial intermediate pole.  They represent the combined effects
   of frame bias, precession and nutation.

3) The fundamental arguments used are as adopted in IERS Conventions
   (2003) and are from Simon et al. (1994) and Souchay et al.  (1999).

4) This is an alternative to the angles-based method, via the ERFA
   function eraFw2xy and as used in eraXys06a for example.  The two
   methods agree at the 1 microarcsecond level (at present), a
   negligible amount compared with the intrinsic accuracy of the
   models.  However, it would be unwise to mix the two methods
   (angles-based and series-based) in a single application.

# References

Capitaine, N., Wallace, P.T. & Chapront, J., 2003, Astron.Astrophys.,
412, 567

Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003), IERS
Technical Note No. 32, BKG

Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M., Francou,
G. & Laskar, J., Astron.Astrophys., 1994, 282, 663

Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M., 1999,
Astron.Astrophys.Supp.Ser. 135, 111

Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
"""
function xy06(day1::AbstractFloat, day2::AbstractFloat)
    Δt = ((day1 - JD2000) + day2)/(100*DAYPERYEAR)

    #  Lunar, solar, and planetary longitudes
    ϕ = deg2rad.(rem.([
        Polynomial(l0_2003A...)(Δt),
        Polynomial(l1_2003A...)(Δt),
        Polynomial( F_2003A...)(Δt),
        Polynomial( D_2003A...)(Δt),
        Polynomial( Ω_2003A...)(Δt)], ARCSECPER2PI)/3600.0)
    append!(ϕ, [
        Polynomial( lme_2003...)(Δt),
        Polynomial( lve_2003...)(Δt),
        Polynomial( lea_2003...)(Δt),
        Polynomial( lma_2003...)(Δt),
        Polynomial( lju_2003...)(Δt),
        Polynomial( lsa_2003...)(Δt),
        Polynomial( lur_2003...)(Δt),
        Polynomial( lne_2003...)(Δt),
        Polynomial( lge_2003...)(Δt)])
    
    #  Polynomial part of precession-nutation
    xypr = SVector(sum(Polynomial(cip_x_2006...)(Δt)),
            sum(Polynomial(cip_y_2006...)(Δt)))

    # !!! The following code can be improved by rearranging the data arrays.

    jaxy = SVector(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1) .+ 1
    jasc = SVector(0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0) .+ 1
    japt = SVector(0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4)

    #  Nutation periodic terms, planetary
    xypl = MVector(0.0, 0.0)
    ialast = length(cip_amplitude_2006)
    for ifreq = length(cip_planetary_2006):-1:1
        sc = sincos(sum(cip_planetary_2006[ifreq].*ϕ))
        ia = cip_pointer_2006[ifreq + length(cip_lunisolar_2006)]
        for i = ialast:-1:ia
            xypl[jaxy[i-ia+1]] += cip_amplitude_2006[i] * sc[jasc[i-ia+1]] * Δt^japt[i-ia+1]
        end
        ialast = ia-1
    end

    #  Nutation periodic terms, luni-solar
    xyls = MVector(0.0, 0.0)
    for ifreq = length(cip_lunisolar_2006):-1:1
        sc = sincos(sum(cip_lunisolar_2006[ifreq].*ϕ[1:5]))
        ia = cip_pointer_2006[ifreq]
        for i = ialast:-1:ia
            xyls[jaxy[i-ia+1]] += cip_amplitude_2006[i] * sc[jasc[i-ia+1]] * Δt^japt[i-ia+1]
        end
        ialast = ia-1
    end

    x, y = deg2rad.((xypr .+ (xyls .+ xypl)./1e6)/3600.0)
    (x = x, y = y)
end

"""
    xys00a(day1::AbstractFloat, day2::AbstractFloat)

For a given TT date, compute the X,Y coordinates of the Celestial
Intermediate Pole and the CIO locator s, using the IAU 2000A
precession-nutation model.

# Input

 - `day1`  -- TT as a 2-part Julian Date (Note 1)
 - `day2`  -- TT as a 2-part Julian Date (Note 1)

# Output

 - `x, y`  -- Celestial Intermediate Pole (Note 2)
 - `s`     -- the CIO locator s (Note 3)

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
   good compromises between resolution and convenience.

2) The Celestial Intermediate Pole coordinates are the x,y components
   of the unit vector in the Geocentric Celestial Reference System.

3) The CIO locator s (in radians) positions the Celestial Intermediate
   Origin on the equator of the CIP.

4) A faster, but slightly less accurate result (about 1 mas for X,Y),
   can be obtained by using instead the eraXys00b function.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function xys00a(day1::AbstractFloat, day2::AbstractFloat)
    #  Form bias-precession-nutation matrix (IAU 2000A) and extract x, y.
    x, y = bpn2xy(pnm00a(day1, day2))
    #  Obtain s
    (x = x, y = y, s = s00(day1, day2, x, y))
end

"""
    xys00b(day1::AbstractFloat, day2::AbstractFloat)

For a given TT date, compute the X,Y coordinates of the Celestial
Intermediate Pole and the CIO locator s, using the IAU 2000B
precession-nutation model.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `x, y`  -- Celestial Intermediate Pole (Note 2)
 - `s`     -- the CIO locator s (Note 3)

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
   good compromises between resolution and convenience.

2) The Celestial Intermediate Pole coordinates are the x,y components
   of the unit vector in the Geocentric Celestial Reference System.

3) The CIO locator s (in radians) positions the Celestial Intermediate
   Origin on the equator of the CIP.

4) The present function is faster, but slightly less accurate (about 1
   mas in X,Y), than the eraXys00a function.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)
"""
function xys00b(day1::AbstractFloat, day2::AbstractFloat)
    #  Form bias-precession-nutation matrix (IAU 2000A) and extract x, y.
    x, y = bpn2xy(pnm00b(day1, day2))
    #  Obtain s
    (x = x, y = y, s = s00(day1, day2, x, y))
end

"""
    xys06a(day1::AbstractFloat, day2::AbstractFloat)

For a given TT date, compute the X,Y coordinates of the Celestial
Intermediate Pole and the CIO locator s, using the IAU 2006 precession
and IAU 2000A nutation models.

# Input

 - `day1`  -- TT as Julian Date (Note 1)
 - `day2`  -- ... as Julian Date

# Output

 - `x, y`  -- Celestial Intermediate Pole (Note 2)
 - `s`     -- the CIO locator s (Note 3)

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
   good compromises between resolution and convenience.

2) The Celestial Intermediate Pole coordinates are the x,y components
   of the unit vector in the Geocentric Celestial Reference System.

3) The CIO locator s (in radians) positions the Celestial Intermediate
   Origin on the equator of the CIP.

4) Series-based solutions for generating X and Y are also available:
   see Capitaine & Wallace (2006) and eraXy06.

# References

Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855

Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
"""
function xys06a(day1::AbstractFloat, day2::AbstractFloat)
    #  Form bias-precession-nutation matrix (IAU 2000A) and extract x, y.
    x, y = bpn2xy(pnm06a(day1, day2))
    #  Obtain s
    (x = x, y = y, s = s06(day1, day2, x, y))
end
