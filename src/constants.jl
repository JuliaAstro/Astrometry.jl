"""
    R2D

Radians to degrees
"""
const R2D = 180/π

"""
    D2R

Degrees to radians
"""
const D2R = π/180

"""
    ASECPERRAD

Arcseconds per radian 

For SOFA name, use [`R2AS`](@ref)
"""
const ASECPERRAD = 180*3600/π

"""
    R2AS

Radians to arcseconds

SOFA name for [`ASECPERRAD`](@ref)
"""
const R2AS = ASECPERRAD

"""
    AS2R

Arcseconds to radians
"""
const AS2R = D2R/3600 # 4.848136811095359935899141e-6

"""
    S2R

Seconds of time to radians
"""
const S2R = 2π/86400 # 7.272205216643039903848712e-5)

"""
    ARCSECPER2PI

Arcseconds per 2pi radians

For SOFA name, use [`TURNAS`](@ref)
"""
const ARCSECPER2PI = 360*3600

"""
    TURNAS

Arcseconds in a full circle

SOFA name for [`ARCSECPER2PI`](@ref)
"""
const TURNAS = ARCSECPER2PI

"""
    MAS2R

Milliarcseconds to radians
"""
const MAS2R = AS2R/1000

"""
    TY

Length of tropical year B1900 (days)
"""
const TY = 365.242198781

"""
    DAYSEC

Seconds per day.
"""
const DAYSEC = 86400

"""
    DAYPERYEAR

Days per Julian year

For SOFA name, use [`JY`](@ref)
"""
const DAYPERYEAR = 365.25

"""
    JY

Days per Julian year

SOFA name for [`DAYPERYEAR`](@ref)
"""
const JY = DAYPERYEAR

"""
    JC

Days per Julian century
"""
const JC = 36525

"""
    JM

Days per Julian millennium
"""
const JM = 365250

"""
    BD1900

J2000.0 - B1900.0 (2415019.81352) in days
"""
const BD1900 = 36524.68648 # Not exposed by SOFA

"""
    JDMIN

Julian day minimum
"""
const JDMIN = -68569.5 # Not exposed by SOFA

"""
    JDMAX

Julian day maximum
"""
const JDMAX = 1e9 # Not exposed by SOFA

"""
    JD2000

Julian day for J2000.0

For SOFA name use [`J00`](@ref)
"""
const JD2000 = 2451545.0

"""
    J00

Julian day for J2000.0

SOFA name for [`JD2000`](@ref)
"""
const J00 = JD2000

"""
    MJD0

Julian day of modified Julian day 0

For SOFA name use [`JM0`](@ref)
"""
const MJD0 = 2400000.5 # day

"""
    JM0

Julian day of modified Julian day 0

SOFA name for [`MJD0`](@ref)
"""
const JM0 = MJD0

"""
    MJD00

Modified Julian day of J2000.0

For SOFA name, use [`JM00`](@ref)
"""
const MJD00 = 51544.5

"""
    JM00

Modified Julian day of J2000.0

SOFA name for [`MJD00`](@ref)
"""
const JM00 = MJD00

"""
    MJD77

Modified Julian day for 1977/01/01

For SOFA name, use [`JM77`](@ref)
"""
const MJD77 = 43144.0

"""
    JM77

Modified Julian day for 1977/01/01

SOFA name for [`MJD77`](@ref)
"""
const JM77 = MJD77

"""
    TT_MINUS_TAI

TT - TAI (s)

For SOFA name use [`TTMTAI`](@ref)
"""
const TT_MINUS_TAI = 32.184 # second

"""
    TTMTAI

TT minus TAI (s)

SOFA name for [`TT_MINUS_TAI`](@ref)
"""
const TTMTAI = TT_MINUS_TAI

"""
    AU

Astronomical unit (m, IAU 2012)
"""
const AU = 1.495978707e11

"""
    CMPS

Speed of light (m/s)
"""
const CMPS = 299792458.0

"""
    AULT

Light time for 1 au (s)
"""
const AULT = AU/CMPS

"""
    C

Speed of light (au per day)
"""
const C = DAYSEC/AULT

"""
    ELG

`L_G = 1 - d(TT)/d(TCG)` (s)
"""
const ELG    = 6.969290134e-10 # second

"""
    TDB0

TDB at TAI 1977/01/01 (s)
"""
const TDB0   = -6.55e-5 # second

"""
    ELB

`L_B = 1 - d(TDB)/d(TCB)` (s)
"""
const ELB    = 1.550519768e-8

"""
    GK

Gaussian gravitational constant

See [`https://en.wikipedia.org/wiki/Gaussian_gravitational_constant`]()
"""
const GK = 0.017202098950 # Not exposed by SOFA

"""
    SCHWARZRADIUS

Schwarzschild radius of the Sun (au)
`2 * 1.32712440041e20 / (2.99792458e8)^2 / 1.49597870700e11`

For SOFA name use [`SRS`](@ref)
"""
const SCHWARZRADIUS = 1.97412574336e-8 # au

"""
    SRS

Schwarzschild radius of the Sun (au)
`2 * 1.32712440041e20 / (2.99792458e8)^2 / 1.49597870700e11`

SOFA name for [`SCHWARZRADIUS`](@ref)
"""
const SRS = SCHWARZRADIUS

"""
    WGS84

WGS84 reference ellipsoid

See also [`GRS80`](@ref) and [`WGS72`].
"""
const WGS84 = :WGS84

"""
    GRS80

GRS80 reference ellipsoid

See also [`WGS84`](@ref) and [`WGS72`].
"""
const GRS80 = :GRS80

"""
    WGS72

WGS72 reference ellipsoid

See also [`WGS84`](@ref) and [`GRS80`].
"""
const WGS72 = :WGS72
