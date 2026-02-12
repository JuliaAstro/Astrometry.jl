#### Astronomy / Time Scales

""" 
    d2dtf(scale::AbstractString, ndp::Int, day1::AbstractFloat, day2::AbstractFloat)

Format for output a 2-part Julian Date (or in the case of UTC a
quasi-JD form that includes special provision for leap seconds).

# Input

 - `scale` -- time scale ID (Note 1)
 - `ndp`   -- resolution (Note 2)
 - `day1`  -- time as a 2-part Julian Date (Notes 3, 4)
 - `day2`  -- time as a 2-part Julian Date (Notes 3, 4)

# Output

 - `date`  -- year, month, day in Gregorian calendar (Note 1, 5)

# Note

1) scale identifies the time scale.  Only the value "UTC" (in upper
   case) is significant, and enables handling of leap seconds (see
   Note 4).

2) ndp is the number of decimal places in the seconds field, and can
   have negative as well as positive values, such as:

   ndp         resolution
   -4            1 00 00
   -3            0 10 00
   -2            0 01 00
   -1            0 00 10
    0            0 00 01
    1            0 00 00.1
    2            0 00 00.01
    3            0 00 00.001

   The limits are platform dependent, but a safe range is -5 to +9.

3) day1+day2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where day1 is the Julian Day Number
   and day2 is the fraction of a day.  In the case of UTC, where the
   use of JD is problematical, special conventions apply: see the next
   note.

4) JD cannot unambiguously represent UTC during a leap second unless
   special measures are taken.  The ERFA internal convention is that
   the quasi-JD day represents UTC days whether the length is 86399,
   86400 or 86401 SI seconds.  In the 1960-1972 era there were smaller
   jumps (in either direction) each time the linear UTC(TAI)
   expression was changed, and these "mini-leaps" are also included in
   the ERFA convention.

5) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale or that are too far in the future to
   be trusted.  See eraDat for further details.

6) For calendar conventions and limitations, see cal2jd.
"""
function d2dtf(scale::AbstractString, ndp::Int, day1::AbstractFloat, day2::AbstractFloat)
    leap = false

    # Provisional calendar date.
    year, month, day, subday = values(jd2cal(day1, day2))

    # If scale is UTC, check for leap second day
    if scale == "UTC"
        # TAI-UTC at  0h on date
        Δt00h = dat(year, month, day, 0.0)

        # TAI-UTC at 12h on date
        Δt12h = dat(year, month, day, 0.5)

        # TAI-UTC at  0h on next date (to detect jumps)
        Δt24h = dat(values(jd2cal(day1+1.5, day2-subday))[1:3]..., 0.0)

        # Check for sudden change in TAI-UTC (seconds)
        Δt = Δt24h - (2*Δt12h - Δt00h)
        
        # if leap second day, scale the fraction of a day into SI.
        leap = abs(Δt) > 0.5
        if leap subday += subday*Δt/SECPERDAY end
    end

    sign, hour, minute, second, subsec = d2tf(ndp, subday)
    
    # Check for rounded time >24 hours
    if hour > 23
        # Check for leap second day
        if leap
            # Check that leap second has passed.
            date = second > 0 ?
                # Use next day but allow for the leap second
                (values(jd2cal(day1+1.5, day2-subday))[1:3]..., 0, 0, 0, subsec) :
                # Use 23 50 60 for the time.
                (year, month, day, 23, 59, 60, subsec)
            if ndp < 0 && second == 60
                date = (values(jd2cal(day1+1.5, day2-subday))[1:3]..., 0, 0, 0, subsec)
            end
        else
            # Use 0h next day
            date = (values(jd2cal(day1+1.5, day2-subday))[1:3]..., 0, 0, 0, subsec)
        end
    else
        date = (year, month, day, hour, minute, second, subsec)
    end
    NamedTuple{(:year, :month, :day, :hour, :minute, :second, :fraction)}(date)
end

"""
    dat(year::Integer, month::Integer, day::Integer, subday::AbstractFloat)

For a given UTC date, calculate Δ(AT) = TAI-UTC.

# Input

 - `year`  -- year (Notes 1 and 2) (in UTC)
 - `month` -- month (Note 2)
 - `day`   -- day (Notes 2 and 3)
 - `fraction` -- fraction of day (Note 4)

# Output

 - `Δt`    -- TAI - UTC (in seconds)

# Note

1) UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper to
   call the function with an earlier date.  If this is attempted, zero
   is returned together with a warning status.

   Because leap seconds cannot, in principle, be predicted in advance,
   a reliable check for dates beyond the valid range is impossible.
   To guard against gross errors, a year five or more after the
   release year of the present function (see the constant IYV) is
   considered dubious.  In this case a warning status is returned but
   the result is computed in the normal way.

   For both too-early and too-late years, the warning status is +1.
   This is distinct from the error status -1, which signifies a year
   so early that JD could not be computed.

2) If the specified date is for a day which ends with a leap second,
   the TAI-UTC value returned is for the period leading up to the leap
   second.  If the date is for a day which begins as a leap second
   ends, the TAI-UTC returned is for the period following the leap
   second.

3) The day number must be in the normal calendar range, for example 1
   through 30 for April.  The "almanac" convention of allowing such
   dates as January 0 and December 32 is not supported in this
   function, in order to avoid confusion near leap seconds.

4) The fraction of day is used only for dates before the introduction
   of leap seconds, the first of which occurred at the end of 1971.
   It is tested for validity (0 to 1 is the valid range) even if not
   used; if invalid, zero is used and status -4 is returned.  For many
   applications, setting fd to zero is acceptable; the resulting error
   is always less than 3 ms (and occurs only pre-1972).

5) The status value returned in the case where there are multiple
   errors refers to the first error detected.  For example, if the
   month and day are 13 and 32 respectively, status -2 (bad month)
   will be returned.  The "internal error" status refers to a case
   that is impossible but causes some compilers to issue a warning.

6) In cases where a valid result is not available, zero is returned.

# References

1) For dates from 1961 January 1 onwards, the expressions from the
   file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.

2) The 5ms timestep at 1961 January 1 is taken from 2.58.1 (p87) of
   the 1992 Explanatory Supplement.
"""
function dat(year::Integer, month::Integer, day::Integer, subday::AbstractFloat)

    @assert 0.0 <= subday <= 1.0 "Fractional day out of range [0-1]."
    @assert year >= DRIFTSECOND[1].year "UTC date is out of range [$(DRIFTSECOND[1].year)-present]."
    if (year > IYV + 5) @warn "UTC date $year-$month-$day is suspect." end

    Δt = typeof(subday)(0.0)
    if year < 1972
        # Find drift offset
        for drift in reverse(DRIFTSECOND)
            if (12*year + month) >= (12*drift.year + drift.month)
                Δt = drift.offset +
                    (cal2jd(year, month, day)[:mjd] - drift.mjd + subday)*drift.rate
                break
            end
        end                     
    else
        # Find leap second
        for leap in reverse(LEAPSECOND)
            if (12*year + month) >= (12*leap.year + leap.month)
                Δt = leap.second
                break
            end
        end
    end
    Δt
end

"""
    dtdb(day1::AbstractFloat, day2::AbstractFloat, ut1::AbstractFloat, eastlon::AbstractFloat,
         u::AbstractFloat, v::AbstractFloat)

An approximation to TDB-TT, the difference between barycentric
dynamical time and terrestrial time, for an observer on the Earth.

The different time scales - proper, coordinate and realized - are
related to each other:

          TAI             <-  physically realized
           :
        offset            <-  observed (nominally +32.184s)
           :
          TT              <-  terrestrial time
           :
  rate adjustment (L_G)   <-  definition of TT
           :
          TCG             <-  time scale for GCRS
           :
    "periodic" terms      <-  eraDtdb  is an implementation
           :
  rate adjustment (L_C)   <-  function of solar-system ephemeris
           :
          TCB             <-  time scale for BCRS
           :
  rate adjustment (-L_B)  <-  definition of TDB
           :
          TDB             <-  TCB scaled to track TT
           :
    "periodic" terms      <-  -eraDtdb is an approximation
           :
          TT              <-  terrestrial time

Adopted values for the various constants can be found in the IERS
Conventions (McCarthy & Petit 2003).

# Input

 - `day1`  -- TBD as part 1 of date (Notes 1-3)
 - `day2`  -- TDB as part 2 of date
 - `ut`    -- universal time (UT1, fraction of one day)
 - `elong` -- longitude (east positive, radians)
 - `u`     -- distance from Earth spin axis (km)
 - `v`     -- distance north of equatorial plane (km)

# Output

 - `tdbtt` -- TDB - TT

# Note

1) The date date1+date2 is a Julian Date, apportioned in any
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

   Although the date is, formally, barycentric dynamical time (TDB),
   the terrestrial dynamical time (TT) can be used with no practical
   effect on the accuracy of the prediction.

2) TT can be regarded as a coordinate time that is realized as an
   offset of 32.184s from International Atomic Time, TAI.  TT is a
   specific linear transformation of geocentric coordinate time TCG,
   which is the time scale for the Geocentric Celestial Reference
   System, GCRS.

3) TDB is a coordinate time, and is a specific linear transformation
   of barycentric coordinate time TCB, which is the time scale for the
   Barycentric Celestial Reference System, BCRS.

4) The difference TCG-TCB depends on the masses and positions of the
   bodies of the solar system and the velocity of the Earth.  It is
   dominated by a rate difference, the residual being of a periodic
   character.  The latter, which is modeled by the present function,
   comprises a main (annual) sinusoidal term of amplitude
   approximately 0.00166 seconds, plus planetary terms up to about 20
   microseconds, and lunar and diurnal terms up to 2 microseconds.
   These effects come from the changing transverse Doppler effect and
   gravitational red-shift as the observer (on the Earth's surface)
   experiences variations in speed (with respect to the BCRS) and
   gravitational potential.

5) TDB can be regarded as the same as TCB but with a rate adjustment
   to keep it close to TT, which is convenient for many applications.
   The history of successive attempts to define TDB is set out in
   Resolution 3 adopted by the IAU General Assembly in 2006, which
   defines a fixed TDB(TCB) transformation that is consistent with
   contemporary solar-system ephemerides.  Future ephemerides will
   imply slightly changed transformations between TCG and TCB, which
   could introduce a linear drift between TDB and TT; however, any
   such drift is unlikely to exceed 1 nanosecond per century.

6) The geocentric TDB-TT model used in the present function is that of
   Fairhead & Bretagnon (1990), in its full form.  It was originally
   supplied by Fairhead (private communications with P.T.Wallace,
   1990) as a Fortran subroutine.  The present C function contains an
   adaptation of the Fairhead code.  The numerical results are
   essentially unaffected by the changes, the differences with respect
   to the Fairhead & Bretagnon original being at the 1e-20 s level.

   The topocentric part of the model is from Moyer (1981) and Murray
   (1983), with fundamental arguments adapted from Simon et al. 1994.
   It is an approximation to the expression ( v / c ) . ( r / c ),
   where v is the barycentric velocity of the Earth, r is the
   geocentric position of the observer and c is the speed of light.

   By supplying zeroes for u and v, the topocentric part of the model
   can be nullified, and the function will return the Fairhead &
   Bretagnon result alone.

7) During the interval 1950-2050, the absolute accuracy is better than
   +/- 3 nanoseconds relative to time ephemerides obtained by direct
   numerical integrations based on the JPL DE405 solar system
   ephemeris.

8) It must be stressed that the present function is merely a model,
   and that numerical integration of solar-system ephemerides is the
   definitive method for predicting the relationship between TCG and
   TCB and hence between TT and TDB.

# References

Fairhead, L., & Bretagnon, P., Astron.Astrophys., 229, 240-247 (1990).

IAU 2006 Resolution 3.

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Moyer, T.D., Cel.Mech., 23, 33 (1981).

Murray, C.A., Vectorial Astrometry, Adam Hilger (1983).

Seidelmann, P.K. et al., Explanatory Supplement to the Astronomical
Almanac, Chapter 2, University Science Books (1992).

Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M., Francou,
G. & Laskar, J., Astron.Astrophys., 282, 663-683 (1994).
"""
function dtdb(day1::AbstractFloat, day2::AbstractFloat, ut1::AbstractFloat, eastlon::AbstractFloat,
              u::AbstractFloat, v::AbstractFloat)
    #  Time since J2000.0 in Julian millenia.
    Δt = ((day1 - JD2000) + day2)/(1000*DAYPERYEAR)

    #  Topocentric terms
    tsol = 2pi*rem(ut1, 1) + eastlon

    #  Fundamental arguments (Simon et al. 1994)
    ϵsun  = deg2rad(rem(Polynomial(ϵsun_1994...)(Δt/3600), 360))
    ϵmsun = deg2rad(rem(Polynomial(ϵmsun_1994...)(Δt/3600), 360))
    D     = deg2rad(rem(Polynomial(D_1994...)(Δt/3600), 360))
    ϵju   = deg2rad(rem(Polynomial(ϵju_1994...)(Δt/3600), 360))
    ϵsa   = deg2rad(rem(Polynomial(ϵsa_1994...)(Δt/3600), 360))

    #  Topocentric terms: Moyer 1981 and Murray 1983
    wt = sum(topo_1994 .* SVector(
        u * sin(tsol + ϵsun - ϵsa),
        u * sin(tsol - 2*ϵmsun),
        u * sin(tsol - D),
        u * sin(tsol + ϵsun - ϵju),
        u * sin(tsol + 2*ϵsun + ϵmsun),
        v * cos(ϵsun + ϵmsun),
        u * sin(tsol - ϵmsun),
        u * sin(tsol + 2*ϵsun),
        v * cos(ϵsun),
        u * sin(tsol)))

    #  Fairhead et al. model
    
    wf = Polynomial(SVector(
        sum(tdb_tt_2003_0[1,:] .*
            sin.(tdb_tt_2003_0[3,:] .+ tdb_tt_2003_0[2,:].*Δt)),
        sum(tdb_tt_2003_1[1,:] .*
            sin.(tdb_tt_2003_1[3,:] .+ tdb_tt_2003_1[2,:].*Δt)),
        sum(tdb_tt_2003_2[1,:] .*
            sin.(tdb_tt_2003_2[3,:] .+ tdb_tt_2003_2[2,:].*Δt)),
        sum(tdb_tt_2003_3[1,:] .*
            sin.(tdb_tt_2003_3[3,:] .+ tdb_tt_2003_3[2,:].*Δt)),
        sum(tdb_tt_2003_4[1,:] .*
            sin.(tdb_tt_2003_4[3,:] .+ tdb_tt_2003_4[2,:].*Δt)))...)(Δt)
    
    wj = sum(mass_plan_1994_0[1,:] .*
             sin.(mass_plan_1994_0[3,:] .+ mass_plan_1994_0[2,:].*Δt)) +
             mass_plan_1994_2 * Δt^2

    wt + wf + wj
end

"""
    dtf2d(scale::AbstractString, year::Int, month::Int, day::Int, hour::Int,
          minute::Int, second::AbstractFloat)

Encode date and time fields into 2-part Julian Date (or in the case of
UTC a quasi-JD form that includes special provision for leap seconds).

# Input

 - `scale`  -- time scale ID (Note 1)
 - `year`   -- year in Gregorian calendar (Note 2)
 - `month`  -- month in Gregorian calendar (Note 2)
 - `day`    -- day in Gregorian calendar (Note 2)
 - `hour`   -- hour
 - `minute` -- minute
 - `second` -- seconds

# Output

 - `day1`   -- part 1 of Julian Date (Notes 3,4)
 -` day2`   -- part 2 of Julian Date (Notes 3,4)

# Note

1) scale identifies the time scale.  Only the value "UTC" (in upper
   case) is significant, and enables handling of leap seconds (see
   Note 4).

2) For calendar conventions and limitations, see eraCal2jd.

3) The sum of the results, d1+d2, is Julian Date, where normally d1 is
   the Julian Day Number and d2 is the fraction of a day.  In the case
   of UTC, where the use of JD is problematical, special

3) The sum of the results, d1+d2, is Julian Date, where normally d1 is
   the Julian Day Number and d2 is the fraction of a day.  In the case
   of UTC, where the use of JD is problematical, special conventions
   apply: see the next note.

4) JD cannot unambiguously represent UTC during a leap second unless
   special measures are taken.  The ERFA internal convention is that
   the quasi-JD day represents UTC days whether the length is 86399,
   86400 or 86401 SI seconds.  In the 1960-1972 era there were smaller
   jumps (in either direction) each time the linear UTC(TAI)
   expression was changed, and these "mini-leaps" are also included in
   the ERFA convention.

5) The warning status "time is after end of day" usually means that
   the sec argument is greater than 60.0.  However, in a day ending in
   a leap second the limit changes to 61.0 (or 59.0 in the case of a
   negative leap second).

6) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale or that are too far in the future to
   be trusted.  See eraDat for further details.

7) Only in the case of continuous and regular time scales (TAI, TT,
   TCG, TCB and TDB) is the result d1+d2 a Julian Date, strictly
   speaking.  In the other cases (UT1 and UTC) the result must be used
   with circumspection; in particular the difference between two such
   results cannot be interpreted as a precise time interval.
"""
function dtf2d(scale::AbstractString, year::Int, month::Int, day::Int, hour::Int,
               minute::Int, second::AbstractFloat)
    # Today's Julian Day number
    julday = sum(cal2jd(year, month, day))

    if scale == "UTC"
        # TAI-UTC at 00h today
        Δt00 = dat(year, month, day, 0.0)

        # TAI-UTC at 12h today
        Δt12 = dat(year, month, day, 0.5)

        # TAI-UTC at 00h next day (to detect jumps)
        Δt24 = dat(values(jd2cal(julday, 1.5))[1:3]..., 0.0)

        # Any sudden change in TAI-UTC between today and tomorrow
        Δt = Δt24 - (2*Δt12 - Δt00)

        # if leap second day, correct the day an final minute lengths
        seclim = hour == 23 && minute == 59 ? 60.0 + Δt : 60.0
    end

    # Validate the time
    @assert 0 <= hour <= 23 "Hour is out of range [0-23]"
    @assert 0 <= minute <= 59 "Minute is out of range [0-59]"
    @assert 0 <= second <= seclim "Second is out of range [0-$seclim]"

    # The time in days
    subday = ((60.0*(60*hour + minute)) + second)/(SECPERDAY + Δt)
    
    (day = julday, fraction = subday)
end

"""
    taitt(day1::AbstractFloat, day2::AbstractFloat)

Time scale transformation: International Atomic Time, TAI, to
Terrestrial Time, TT.

# Input

 - `day1`  -- TAI as part 1 of Julian Date
 - `day2`  -- TAI as part 2 of Julian Date

# Output

 - `tt1`   -- TT as part 1 of Julian Date
 - `tt2`   -- TT as part 2 of Julian Date

# Note

1) day1+day2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where day1 is the Julian Day Number
   and day2 is the fraction of a day.  The returned (day1,day2) follow
   suit.

# References

1) McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
   Technical Note No. 32, BKG (2004)

2) Explanatory Supplement to the Astronomical Almanac, P. Kenneth
   Seidelmann (ed), University Science Books (1992)
"""
function taitt(day1::AbstractFloat, day2::AbstractFloat)
    abs(day1) > abs(day2) ? (day = day1, fraction = day2 + TT_MINUS_TAI/SECPERDAY) :
        (day = day1 + TT_MINUS_TAI/SECPERDAY, fraction = day2)
end

"""
    taiut1(day1::AbstractFloat, day2::AbstractFloat, Δt::AbstractFloat)

Time scale transformation: International Atomic Time, TAI, to
Universal Time, UT1.

# Input

 - `day1`  -- TAI as part 1 of Julian Date
 - `day2`  -- TAI as part 2 of Julian Date
 - `dta`   -- UT1-TAI in seconds

# Output

 - `ut11`  -- UT1 as part 1 of Julian Date
 - `ut12`  -- UT1 as part 2 of Julian Date

# Note

1) day1+day2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where tai1 is the Julian Day Number
   and tai2 is the fraction of a day.  The returned (day1,day2) follow
   suit.

2) The argument dta, i.e. UT1-TAI, is an observed quantity, and is
   available from IERS tabulations.

# References

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992)
"""
function taiut1(day1::AbstractFloat, day2::AbstractFloat, Δt::AbstractFloat)
    abs(day1) > abs(day2) ? (day = day1, fraction = day2 + Δt/SECPERDAY) :
        (day = day1 + Δt/SECPERDAY, fraction = day2)
end

"""
    taiutc(day1::AbstractFloat, day2::AbstractFloat)

Time scale transformation: International Atomic Time, TAI, to
Coordinated Universal Time, UTC.

# Input

 - `day1`  -- TAI as part 1 of Julian Date (Note 1)
 - `day2`  -- TAI as part 2 of Julian Date

# Output

 - `utc1`  -- UTC as part 1 of quasi Julian Date (Notes 1-3)
 - `utc2`  -- UTC as part 2 of quasi Julian Date

# Note

1) day1+day2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where day1 is the Julian Day Number
   and day2 is the fraction of a day.  The returned (utc1, utc2) form
   an analogous pair, except that a special convention is used, to
   deal with the problem of leap seconds - see the next note.

2) JD cannot unambiguously represent UTC during a leap second unless
   special measures are taken.  The convention in the present function
   is that the JD day represents UTC days whether the length is 86399,
   86400 or 86401 SI seconds.  In the 1960-1972 era there were smaller
   jumps (in either direction) each time the linear UTC(TAI)
   expression was changed, and these "mini-leaps" are also included in
   the ERFA convention.

3) The function eraD2dtf can be used to transform the UTC quasi-JD
   into calendar date and clock time, including UTC leap second
   handling.

4) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale or that are too far in the future to
   be trusted.  See eraDat for further details.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992)
"""
function taiutc(day1::AbstractFloat, day2::AbstractFloat)
    utc1, utc2 = tai1, tai2 = abs(day1) >= abs(day2) ? (day1, day2) : (day2, day1)

    for j=1:3
        #  Guess UTC to TAI and adjust guessed UTC
        utc2 += sum((tai1, tai2) .- values(utctai(utc1, utc2)))
    end
    abs(day1) >= abs(day2) ? (day = utc1, fraction = utc2) : (day = utc2, fraction = utc1)
end

"""
    tcbtdb(day1::AbstractFloat, day2::AbstractFloat)

Time scale transformation: Barycentric Coordinate Time, TCB, to
Barycentric Dynamical Time, TDB.

# Input

 - `day1`  -- TCB as part 1 of Julian Date
 - `day2`  -- TCB as part 2 of Julian Date

# Output

 - `tdb1`  -- TDB as part 1 of Julian Date
 - `tdb2`  -- TDB as part 2 of Julian Date

# Note

1) day1+day2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where day1 is the Julian Day Number
   and day2 is the fraction of a day.  The returned (tdb1,tdb2) follow
   suit.

2) The 2006 IAU General Assembly introduced a conventional linear
   transformation between TDB and TCB.  This transformation
   compensates for the drift between TCB and terrestrial time TT, and
   keeps TDB approximately centered on TT.  Because the relationship
   between TT and TCB depends on the adopted solar system ephemeris,
   the degree of alignment between TDB and TT over long intervals will
   vary according to which ephemeris is used.  Former definitions of
   TDB attempted to avoid this problem by stipulating that TDB and TT
   should differ only by periodic effects.  This is a good description
   of the nature of the relationship but eluded precise mathematical
   formulation.  The conventional linear relationship adopted in 2006
   sidestepped these difficulties whilst delivering a TDB that in
   practice was consistent with values before that date.

3) TDB is essentially the same as Teph, the time argument for the JPL
   solar system ephemerides.

# References

IAU 2006 Resolution B3
"""
function tcbtdb(day1::AbstractFloat, day2::AbstractFloat)
    abs(day1) > abs(day2) ?
        (day = day1, fraction = day2 + TDB0/SECPERDAY -
         ELB*((day1 - MJD0 - MJD77) + (day2 - TT_MINUS_TAI/SECPERDAY))) :
         (day = day1 + TDB0/SECPERDAY -
          ELB*((day1 - MJD0 - MJD77) + (day1 - TT_MINUS_TAI/SECPERDAY)),
          fraction = day2)
end

"""
    tcgtt(day1::AbstractFloat, day2::AbstractFloat)

Time scale transformation: Geocentric Coordinate Time, TCG, to
Terrestrial Time, TT.

# Input

 - `day1`  -- TCG as part 1 of Julian Date
 - `day2`  -- TCG as part 2 of Julian Date

# Output

 - `tt1`   -- TT as part 1 of Julian Date
 - `tt1`   -- TT as part 2 of Julian Date

# Note

1) day1+day2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where day1 is the Julian Day Number
   and day2 is the fraction of a day.  The returned (tt1, tt2) follow
   suit.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

IAU 2000 Resolution B1.9
"""
function tcgtt(day1::AbstractFloat, day2::AbstractFloat)
    abs(day1) > abs(day2) ?
        (day = day1, fraction = day2 -
        ELG*((day1 - MJD0) + (day2 - MJD77 - TT_MINUS_TAI/SECPERDAY))) :
        (day = day1 -
        ELG*((day2 - MJD0) + (day1 - MJD77 - TT_MINUS_TAI/SECPERDAY)),
    fraction = day2)
end

"""
    tdbtcb(day1::AbstractFloat, day2::AbstractFloat)

Time scale transformation: Barycentric Dynamical Time, TDB, to
Barycentric Coordinate Time, TCB.

# Input

 - `day1`   -- TDB as part 1 of Julian Date
 - `day2`   -- TDB as part 2 of Julian Date

# Output

 - `tcb1`   -- TCB as part 1 of Julian Date
 - `tcb2`   -- TCB as part 2 of Julian Date

# Note

1) day1+day2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where day1 is the Julian Day Number
   and day2 is the fraction of a day.  The returned tcb1,tcb2 follow
   suit.

2) The 2006 IAU General Assembly introduced a conventional linear
   transformation between TDB and TCB.  This transformation
   compensates for the drift between TCB and terrestrial time TT, and
   keeps TDB approximately centered on TT.  Because the relationship
   between TT and TCB depends on the adopted solar system ephemeris,
   the degree of alignment between TDB and TT over long intervals will
   vary according to which ephemeris is used.  Former definitions of
   TDB attempted to avoid this problem by stipulating that TDB and TT
   should differ only by periodic effects.  This is a good description
   of the nature of the relationship but eluded precise mathematical
   formulation.  The conventional linear relationship adopted in 2006
   sidestepped these difficulties whilst delivering a TDB that in
   practice was consistent with values before that date.

3) TDB is essentially the same as Teph, the time argument for the JPL
   solar system ephemerides.

# References

IAU 2006 Resolution B3
"""
function tdbtcb(day1::AbstractFloat, day2::AbstractFloat)
    abs(day1) > abs(day2) ?
        (day = day1, fraction = day2 - TDB0/SECPERDAY -
        ELB/(1.0 - ELB)*((MJD0 + MJD77 - day1) -
                         (day2 - (TDB0 + TT_MINUS_TAI)/SECPERDAY))) :
        (day = day1 - TDB0/SECPERDAY -
        ELB/(1.0 - ELB)*((MJD0 + MJD77 - day2) -
                         (day1 - (TDB0 + TT_MINUS_TAI)/SECPERDAY)),
    fraction = day2)
end

"""
    tdbtt(day1::AbstractFloat, day2::AbstractFloat, dtr::AbstractFloat)

Time scale transformation: Barycentric Dynamical Time, TDB, to
Terrestrial Time, TT.

# Input

 - `day1`  -- TDB as part 1 of Julian Date
 - `day2`  -- TDB as part 2 of Julian Date
 - `dtr`   -- TDB-TT in seconds

# Output

 - `tt1`   -- TT as part 1 of Julian Date
 - `tt2`   -- TT as part 2 of Julian Date

# Note

1) day1+day2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where day1 is the Julian Day Number
   and day2 is the fraction of a day.  The returned tt1,tt2 follow
   suit.

2) The argument dtr represents the quasi-periodic component of the GR
   transformation between TT and TCB.  It is dependent upon the
   adopted solar-system ephemeris, and can be obtained by numerical
   integration, by interrogating a precomputed time ephemeris or by
   evaluating a model such as that implemented in the ERFA function
   eraDtdb.  The quantity is dominated by an annual term of 1.7 ms
   amplitude.

3) TDB is essentially the same as Teph, the time argument for the JPL
   solar system ephemerides.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

IAU 2006 Resolution 3
"""
function tdbtt(day1::AbstractFloat, day2::AbstractFloat, dtr::AbstractFloat)
    abs(day1) > abs(day2) ? (day = day1, fraction = day2 - dtr/SECPERDAY) :
        (day = day1 - dtr/SECPERDAY, fraction = day2)
end

"""
    tttai(day1::AbstractFloat, day2::AbstractFloat)

Time scale transformation: Terrestrial Time, TT, to International
Atomic Time, TAI.

# Input

 - `day1`  -- TT as part 1 of Julian Date
 - `day2`  -- TT as part 2 of Julian Date

# Output

 - `tai1`  -- TAI as part 1 of Julian Date
 - `tai2`  -- TAI as part 2 of Julian Date

# Note

1) day1+day2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where day1 is the Julian Day Number
   and day2 is the fraction of a day.  The returned tai1,tai2 follow
   suit.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992)
"""
function tttai(day1::AbstractFloat, day2::AbstractFloat)
    abs(day1) > abs(day2) ?
        (day = day1, frcation = day2 - TT_MINUS_TAI/SECPERDAY) :
        (day = day1 - TT_MINUS_TAI/SECPERDAY, fraction = day2)
end

"""
    tttcg(day1::AbstractFloat, day2::AbstractFloat)

Time scale transformation: Terrestrial Time, TT, to Geocentric
Coordinate Time, TCG.

# Input

 - `day1`  -- TT as part 1 of Julian Date
 - `day2`  -- TT as part 2 of Julian Date

# Output

 - `tcg1`  -- TCG as part 1 of Julian Date
 - `tcg2`  -- TCG as part 2 of Julian Date

# Note

1) day1+day2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where day1 is the Julian Day Number
   and day2 is the fraction of a day.  The returned tcg1,tcg2 follow
   suit.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

IAU 2000 Resolution B1.9
"""
function tttcg(day1::AbstractFloat, day2::AbstractFloat)
    abs(day1) > abs(day2) ?
        (day = day1, fraction = day2 +
        ELG/(1.0 - ELG)*((day1 - MJD0) +
                         (day2 - MJD77 - TT_MINUS_TAI/SECPERDAY))) :
        (day = day1 +
        ELG/(1.0 - ELG)*((day2 - MJD0) +
                         (day1 - MJD77 - TT_MINUS_TAI/SECPERDAY)),
    fraction = day2)
end

"""
    tttdb(day1::AbstractFloat, day2::AbstractFloat, dtr::AbstractFloat)

Time scale transformation: Terrestrial Time, TT, to Barycentric
Dynamical Time, TDB.

# Input

 - `day1`  -- TT as part 1 of Julian Date
 - `day2`  -- TT as part 2 of Julian Date
 - `dtr`   -- TDB-TT in seconds

# Output

 - `tdb1`  -- TDB as part 1 of Julian Date
 - `tdb2`  -- TDB as part 2 of Julian Date

# Note

1) day1+day2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where day1 is the Julian Day Number
   and day2 is the fraction of a day.  The returned tdb1,tdb2 follow
   suit.

2) The argument dtr represents the quasi-periodic component of the GR
   transformation between TT and TCB.  It is dependent upon the
   adopted solar-system ephemeris, and can be obtained by numerical
   integration, by interrogating a precomputed time ephemeris or by
   evaluating a model such as that implemented in the ERFA function
   eraDtdb.  The quantity is dominated by an annual term of 1.7 ms
   amplitude.

3) TDB is essentially the same as Teph, the time argument for the JPL
   solar system ephemerides.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

IAU 2006 Resolution 3
"""
function tttdb(day1::AbstractFloat, day2::AbstractFloat, dtr::AbstractFloat)
    abs(day1) > abs(day2) ?
        (day = day1, fraction = day2 + dtr/SECPERDAY) :
        (day = day1 + dtr/SECPERDAY, fraction = day2)
end

"""
    ttut1(day1::AbstractFloat, day2::AbstractFloat, dt::AbstractFloat)

Time scale transformation: Terrestrial Time, TT, to Universal Time,
UT1.

# Input

 - `day1`  -- TT as part 1 of Julian Date
 - `day2`  -- TT as part 2 of Julian Date
 - `dt`    -- TT-UT1 in seconds

# Output

 - `ut11`  -- UT1 as part 1 of Julian Date
 - `ut12`  -- UT1 as part 2 of Julian Date

# Note

1) day1+day2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where day1 is the Julian Day Number
   and day2 is the fraction of a day.  The returned ut11,ut12 follow
   suit.

2) The argument dt is classical Delta T.

# References

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992)
"""
function ttut1(day1::AbstractFloat, day2::AbstractFloat, dt::AbstractFloat)
    abs(day1) > abs(day2) ?
        (day = day1, fraction = day2 - dt/SECPERDAY) :
        (day = day1 - dt/SECPERDAY, fraction = day2)
end

"""
    ut1tai(day1::AbstractFloat, day2::AbstractFloat, dta::AbstractFloat)

Time scale transformation: Universal Time, UT1, to International
Atomic Time, TAI.

# Input

 - `day1`  -- UT1 as part 1 of Julian Date
 - `day2`  -- UT1 as part 2 of Julian Date
 - `dta`   -- UT1-TAI in seconds

# Output

 - `tai1`  -- TAI as part 1 of Julian Date
 - `tai2`  -- TAI as part 2 of Julian Date

# Note

1) day1+day2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where day1 is the Julian Day Number
   and day2 is the fraction of a day.  The returned tai1,tai2 follow
   suit.

2) The argument dta, i.e. UT1-TAI, is an observed quantity, and is
   available from IERS tabulations.

# References

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992)
"""
function ut1tai(day1::AbstractFloat, day2::AbstractFloat, dta::AbstractFloat)
    abs(day1) > abs(day2) ?
        (day = day1, fraction = day2 - dta/SECPERDAY) :
        (day = day1 - dta/SECPERDAY, fraction = day2)
end

"""
    ut1tt(day1::AbstractFloat, day2::AbstractFloat, dt::AbstractFloat)

Time scale transformation: Universal Time, UT1, to Terrestrial Time,
TT.

# Input

 - `day1`  -- UT1 as part 1 Julian Date (Note 1)
 - `day2`  -- UT1 as part 2 Julian Date
 - `dt`    -- TT-UT1 in seconds (Note 2)

# Output

 - `tt1`   -- TT as part 1 Julian Date
 - `tt2`   -- TT as part 2 Julian Date

# Note

1) day1+day2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where day1 is the Julian Day Number
   and day2 is the fraction of a day.  The returned tt1,tt2 follow
   suit.

2) The argument dt is classical Delta T.

# References

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992)
"""
function ut1tt(day1::AbstractFloat, day2::AbstractFloat, dt::AbstractFloat)
    abs(day1) > abs(day2) ?
        (day = day1, fraction = day2 + dt/SECPERDAY) :
        (day = day1 + dt/SECPERDAY, fraction = day2)
end

"""
    ut1utc(day1::AbstractFloat, day2::AbstractFloat, duts::AbstractFloat)

Time scale transformation: Universal Time, UT1, to Coordinated
Universal Time, UTC.

# Input

 - `ut1`   -- UT1 as part 1 of Julian Date (Note 1)
 - `ut2`   -- UT1 as part 2 of Julian Date
 - `duts`  -- Delta UT1: UT1-UTC in seconds (Note 2)

# Output

 - `utc1`  -- UTC as part 1 of quasi Julian Date (Notes 3,4)
 - `utc2`  -- UTC as part 2 of quasi Julian Date (Notes 3,4)

# Note

1) day1+day2 is Julian Date, apportioned in any convenient way
   between the two arguments, for example where day1 is the Julian
   Day Number and day2 is the fraction of a day.  The returned utc1
   and utc2 form an analogous pair, except that a special convention
   is used, to deal with the problem of leap seconds - see Note 3.

2) Delta UT1 can be obtained from tabulations provided by the
   International Earth Rotation and Reference Systems Service.  The
   value changes abruptly by 1s at a leap second;  however, close to
   a leap second the algorithm used here is tolerant of the "wrong"
   choice of value being made.

3) JD cannot unambiguously represent UTC during a leap second unless
   special measures are taken.  The convention in the present function
   is that the returned quasi-JD UTC1+UTC2 represents UTC days whether
   the length is 86399, 86400 or 86401 SI seconds.

4) The function eraD2dtf can be used to transform the UTC quasi-JD
   into calendar date and clock time, including UTC leap second
   handling.

5) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale or that are too far in the future to
   be trusted.  See eraDat for further details.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992)
"""
function ut1utc(day1::AbstractFloat, day2::AbstractFloat, duts::AbstractFloat)
    #  See if the UT1 can possibly be in a leap-second day
    utc1, utc2 = abs(day1) >= abs(day2) ? (day1, day2) : (day2, day1)
    day1, dats1 = utc1, 0.0

    for i=-1:3
        year, month, day, frac = jd2cal(day1, day2+i)
        dats2 = dat(year, month, day, 0.0)
        if i == -1 dats1 = dats2 end
        ddats = dats2 - dats1
        if abs(ddats) >= 0.5
            #  Yes. Leap second nearby: ensure UT1-UTC is 'before' value.
            if ddats * duts >= 0 duts -= ddats end

            #  UT1 for the start of the UTC day that ends in a leap.
            #  Is th UT! after this time?
            du = sum((utc1, utc2 + 1.0 - duts/SECPERDAY) .- cal2jd(year, month, day))

            if du > 0
                #  Yes. fraction of hte current UTC day that has elapsed.
                frac = du*SECPERDAY/(SECPERDAY + ddats)

                #  Ramp UT1-UTC to bring about ERFA's JD(UTC) convention.
                duts += ddats * (frac <= 1.0 ? frac : 1.0)
            end
            break
        end
        dats1 = dats2
    end

    #  Subtract the (possibly adjusted) UT1-UTC from UT1 to give UTC.
    abs(day1) > abs(day2) ?
        (day = utc1, fraction = utc2 - duts/SECPERDAY) :
        (day = utc1 - duts/SECPERDAY, fraction = utc2)
end

"""
    utctai(day1::AbstractFloat, day2::AbstractFloat)

Time scale transformation: Coordinated Universal Time, UTC, to
International Atomic Time, TAI.

# Input

 - `day1`  -- UTC as part 1 of quasi Julian Date (Notes 1-4)
 - `day2`  -- UTC as part 2 of quasi Julian Date

# Output

 - `tai1`  -- TAI as part 1 of Julian Date (Note 5)
 - `tai2`  -- TAI as part 2 of Julian Date

# Note

1) day1+day2 is quasi Julian Date (see Note 2), apportioned in any
   convenient way between the two arguments, for example where day1
   is the Julian Day Number and day2 is the fraction of a day.

2) JD cannot unambiguously represent UTC during a leap second unless
   special measures are taken.  The convention in the present function
   is that the JD day represents UTC days whether the length is 86399,
   86400 or 86401 SI seconds.  In the 1960-1972 era there were smaller
   jumps (in either direction) each time the linear UTC(TAI)
   expression was changed, and these "mini-leaps" are also included in
   the ERFA convention.

3) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale or that are too far in the future to
   be trusted.  See eraDat for further details.

4) The function eraDtf2d converts from calendar date and time of day
   into 2-part Julian Date, and in the case of UTC implements the
   leap-second-ambiguity convention described above.

5) The returned TAI1,TAI2 are such that their sum is the TAI Julian
   Date.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992)
"""
function utctai(day1::AbstractFloat, day2::AbstractFloat)

    utc1, utc2 = abs(day1) >= abs(day2) ? (day1, day2) : (day2, day1)
    year, month, day, frac = jd2cal(utc1, utc2)
    #  Get TAI-UTC at 0 hours today
    dat00 = dat(year, month, day, 0.0)

    #  Get TAI-UTC at 12 hours today (to detect drift).
    dat12 = dat(year, month, day, 0.5)

    #  Get TAI-UTC at 0 hours tomorrow (to detect jumps).
    year1, month1, day1, frac1 = jd2cal(utc1 + 1.5, utc2 - frac)
    dat24 = dat(year1, month1, day1, 0.0)

    #  Separate TAI-UTC change into per-day (DKOD) and any jump (DLEAP).
    dlod = 2*(dat12 - dat00)
    dleap = dat24 - (dat00 + dlod)

    #  Remove any scaling applied to spread leap into preceding day.
    frac *= 1.0 + dleap/SECPERDAY

    #  Scale from (pre-1972) UTC second to SI seconds.
    frac *= 1.0 + dlod/SECPERDAY
    
    #  Today's calendar date to JD
    today1, today2 = cal2jd(year, month, day)
    
    abs(day1) > abs(day2) ?
        (day = utc1, fraction = today1 - utc1 + today2 + frac + dat00/SECPERDAY) :
        (day = today1 - utc1 + today2 + frac + dat00/SECPERDAY, fraction = utc1)
end

"""
    utcut1(day1::AbstractFloat, day2::AbstractFloat, dut1::AbstractFloat)

Time scale transformation: Coordinated Universal Time, UTC, to
Universal Time, UT1.

# Input

 - `day1`  -- UTC as part 1 of quasi Julian Date (Notes 1-4)
 - `day2`  -- UTC as part 2 of quasi Julian Date
 - `dut1`  -- Delta UT1 = UT1-UTC in seconds (Note 5)

# Output

 - `ut11`  -- UT1 as part 1 of Julian Date
 - `ut12`  -- UT1 as part 2 of Julian Date

# Note

1) day1+day2 is quasi Julian Date (see Note 2), apportioned in any
   convenient way between the two arguments, for example where day1 is
   the Julian Day Number and day2 is the fraction of a day.

2) JD cannot unambiguously represent UTC during a leap second unless
   special measures are taken.  The convention in the present function
   is that the JD day represents UTC days whether the length is 86399,
   86400 or 86401 SI seconds.

3) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale or that are too far in the future to
   be trusted.  See eraDat for further details.

4) The function eraDtf2d converts from calendar date and time of day
   into 2-part Julian Date, and in the case of UTC implements the
   leap-second-ambiguity convention described above.

5) Delta UT1 can be obtained from tabulations provided by the
   International Earth Rotation and Reference Systems Service.  It is
   the caller's responsibility to supply a dut1 argument containing
   the UT1-UTC value that matches the given UTC.

6) The returned ut11,ut12 are such that their sum is the UT1 Julian
   Date.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Explanatory Supplement to the Astronomical Almanac, P. Kenneth
Seidelmann (ed), University Science Books (1992)
"""
function utcut1(day1::AbstractFloat, day2::AbstractFloat, dut1::AbstractFloat)
    
    #  Loop up TAI-UTC.
    year, month, day, frac = jd2cal(day1, day2)

    #  Form UT1-TAI and UTC to TAI to UT1.
    NamedTuple{(:day, :fraction)}(
    taiut1(utctai(day1, day2)..., dut1 - dat(year, month, day, 0.0)))
end
