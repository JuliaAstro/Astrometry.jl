#   IAU 2000 Model Common Constants

#   Geodetic coefficients
const wgs72_radius = 6378135.0
const wgs72_oblate = 1.0 / 298.26

const grs80_radius = 6378137.0
const grs80_oblate = 1.0 / 298.257222101

const wgs84_radius = 6378137.0
const wgs84_oblate = 1.0 / 298.257223563

#   Coefficients of IAU 1982 GMST-UT1 model
const gmst_1982 = (24110.54841-SECPERDAY/2.0, 8640184.812866, 0.093104, -6.2e-6)

#   Equation of the equinoxes (1994)
const equinox_1994 = [0.00264, 0.000063]

#   Longitude of the mean ascending node of the lunar orbit (1994)
const l_1994 = [450160.280, -482890.539, 7.455, 0.008]

#   Initial obliquity of the ecliptic in arcseconds
const ϵ0_2000 = 84381.448

#   Earth rotation angle
const era_2000 = (0.7790572732640, 0.00273781191135448)

#   Greenwich Mean Sidereal Time
const gmst_2000 = (0.014506, 4612.15739966, 1.39667721, -0.00009344, 0.00001882)

#   TIO locator s'
const tio_2000 = -4.7e-5

#   Precession model constants in mas (Lieske et al. 1977)
const ψ_1977 = (    0.0, 5038.7784, -1.07259, -0.001147)
const ω_1977 = (ϵ0_2000,    0.0,     0.05127, -0.007726)
const χ_1977 = (    0.0,   10.5526, -2.38064, -0.001125)

#   Frame bias model constants in mas
const icrs_ra_2000 = -0.0146
const ψ_bias_2000  = -0.041775
const ϵ_bias_2000  = -0.0068192

#   Precession model correction constants in mas
const ψ_corr_2000 = -0.29965
const ϵ_corr_2000 = -0.02524

#   Secular variation correction factor
const j2_corr_2000 = -2.7774e-6

#   P03 adjustments (Wallace & Capitaine, 2006, Eq. 5)
const p03_2000 = 0.4697e-6

#   Matrix element for orienting the analytical model to DE405
const iau_2000_bcrs = SMatrix{3,3}(
     1.0,            -0.000000230286, 0.0           ,
     0.000000211284,  0.917482137087, 0.397776982902,
    -0.000000091603, -0.397776982902, 0.917482137087)

#   Matrix element for galactic to ICRS coordinates
const r_gal_icrs = SMatrix{3,3}(
    -0.054875560416215368492398900454, +0.494109427875583673525222371358, -0.867666149019004701181616534570,
    -0.873437090234885048760383168409, -0.444829629960011178146614061616, -0.198076373431201528180486091412,
    -0.483835015548713226831774175116, +0.746982244497218890527388004556, +0.455983776175066922272100478348)

#   Coefficients for fundamental arguments:
#
#   . Powers of time in Julian centuries
#   . Units are degrees.
#
#   Moon's mean longitude (wrt mean equinox and ecliptic of date)
#   Simon et al. (1994).
const λmoon_1994 = (218.31665436, 481267.88123421, -0.0015786, 1/538841.0, -1/65194000.0)

#   Moon's mean elongation
const dmoon_1998 = (297.8501921,  445267.1114034, -0.0018819, 1/545868.0, 1/113065000.0)

#   Sun's mean anomaly
const lsun_1998 = (357.5291092, 35999.0502909, -0.0001536, 1/24490000.0, 0.0)

#   Moon's mean anomaly
const lmoon_1998 = (134.9633964, 477198.8675055, 0.0087414, 1/69699.0, -1/14712000.0)

#   Mean distance of the Moon from its ascending node
const fmoon_1998 = (93.2720950, 483202.0175233, -0.0036539, 1/3526000.0, 1/863310000.0)

#   Other arguments

#   Meeus A_1, due to Venus (deg)
const a_1 = (119.75, 131.849)

#   Meeus A_2, due to Jupiter (deg)
const a_2 = (53.09, 479264.290)

#   Meeus A_3, due to sidereal motion of the Moon in longitude (deg)
const a_3 = (313.45, 481266.484)

#   Coefficients for Meeus "additive terms" (deg)
const a_l = (0.003958, 0.001962, 0.000318)
const a_b = (-0.002235, 0.000382, 0.000175, 0.000175, 0.000127, -0.000115)

#   Fixed term in distance (m)
const r0 = 385000560.0

#   Coefficients for (dimensionless) E factor
const efac = (1., -0.002516, -0.0000074)

#   Coefficients for 

#   Coefficients for Moon longitude and distance series
const lr_1998 = [
    PeriodicTerms([0,  0,  1,  0], [ 6.288774, -20905355.0]),
    PeriodicTerms([2,  0, -1,  0], [ 1.274027,  -3699111.0]),
    PeriodicTerms([2,  0,  0,  0], [ 0.658314,  -2955968.0]),
    PeriodicTerms([0,  0,  2,  0], [ 0.213618,   -569925.0]),
    PeriodicTerms([0,  1,  0,  0], [-0.185116,     48888.0]),
    PeriodicTerms([0,  0,  0,  2], [-0.114332,     -3149.0]),
    PeriodicTerms([2,  0, -2,  0], [ 0.058793,    246158.0]),
    PeriodicTerms([2, -1, -1,  0], [ 0.057066,   -152138.0]),
    PeriodicTerms([2,  0,  1,  0], [ 0.053322,   -170733.0]),
    PeriodicTerms([2, -1,  0,  0], [ 0.045758,   -204586.0]),
    PeriodicTerms([0,  1, -1,  0], [-0.040923,   -129620.0]),
    PeriodicTerms([1,  0,  0,  0], [-0.034720,    108743.0]),
    PeriodicTerms([0,  1,  1,  0], [-0.030383,    104755.0]),
    PeriodicTerms([2,  0,  0, -2], [ 0.015327,     10321.0]),
    PeriodicTerms([0,  0,  1,  2], [-0.012528,         0.0]),
    PeriodicTerms([0,  0,  1, -2], [ 0.010980,     79661.0]),
    PeriodicTerms([4,  0, -1,  0], [ 0.010675,    -34782.0]),
    PeriodicTerms([0,  0,  3,  0], [ 0.010034,    -23210.0]),
    PeriodicTerms([4,  0, -2,  0], [ 0.008548,    -21636.0]),
    PeriodicTerms([2,  1, -1,  0], [-0.007888,     24208.0]),
    PeriodicTerms([2,  1,  0,  0], [-0.006766,     30824.0]),
    PeriodicTerms([1,  0, -1,  0], [-0.005163,     -8379.0]),
    PeriodicTerms([1,  1,  0,  0], [ 0.004987,    -16675.0]),
    PeriodicTerms([2, -1,  1,  0], [ 0.004036,    -12831.0]),
    PeriodicTerms([2,  0,  2,  0], [ 0.003994,    -10445.0]),
    PeriodicTerms([4,  0,  0,  0], [ 0.003861,    -11650.0]),
    PeriodicTerms([2,  0, -3,  0], [ 0.003665,     14403.0]),
    PeriodicTerms([0,  1, -2,  0], [-0.002689,     -7003.0]),
    PeriodicTerms([2,  0, -1,  2], [-0.002602,         0.0]),
    PeriodicTerms([2, -1, -2,  0], [ 0.002390,     10056.0]),
    PeriodicTerms([1,  0,  1,  0], [-0.002348,      6322.0]),
    PeriodicTerms([2, -2,  0,  0], [ 0.002236,     -9884.0]),
    PeriodicTerms([0,  1,  2,  0], [-0.002120,      5751.0]),
    PeriodicTerms([0,  2,  0,  0], [-0.002069,         0.0]),
    PeriodicTerms([2, -2, -1,  0], [ 0.002048,     -4950.0]),
    PeriodicTerms([2,  0,  1, -2], [-0.001773,      4130.0]),
    PeriodicTerms([2,  0,  0,  2], [-0.001595,         0.0]),
    PeriodicTerms([4, -1, -1,  0], [ 0.001215,     -3958.0]),
    PeriodicTerms([0,  0,  2,  2], [-0.001110,         0.0]),
    PeriodicTerms([3,  0, -1,  0], [-0.000892,      3258.0]),
    PeriodicTerms([2,  1,  1,  0], [-0.000810,      2616.0]),
    PeriodicTerms([4, -1, -2,  0], [ 0.000759,     -1897.0]),
    PeriodicTerms([0,  2, -1,  0], [-0.000713,     -2117.0]),
    PeriodicTerms([2,  2, -1,  0], [-0.000700,      2354.0]),
    PeriodicTerms([2,  1, -2,  0], [ 0.000691,         0.0]),
    PeriodicTerms([2, -1,  0, -2], [ 0.000596,         0.0]),
    PeriodicTerms([4,  0,  1,  0], [ 0.000549,     -1423.0]),
    PeriodicTerms([0,  0,  4,  0], [ 0.000537,     -1117.0]),
    PeriodicTerms([4, -1,  0,  0], [ 0.000520,     -1571.0]),
    PeriodicTerms([1,  0, -2,  0], [-0.000487,     -1739.0]),
    PeriodicTerms([2,  1,  0, -2], [-0.000399,         0.0]),
    PeriodicTerms([0,  0,  2, -2], [-0.000381,     -4421.0]),
    PeriodicTerms([1,  1,  1,  0], [ 0.000351,         0.0]),
    PeriodicTerms([3,  0, -2,  0], [-0.000340,         0.0]),
    PeriodicTerms([4,  0, -3,  0], [ 0.000330,         0.0]),
    PeriodicTerms([2, -1,  2,  0], [ 0.000327,         0.0]),
    PeriodicTerms([0,  2,  1,  0], [-0.000323,      1165.0]),
    PeriodicTerms([1,  1, -1,  0], [ 0.000299,         0.0]),
    PeriodicTerms([2,  0,  3,  0], [ 0.000294,         0.0]),
    PeriodicTerms([2,  0, -1, -2], [ 0.000000,      8752.0])]

#   Coefficients for Moon latitude series
const b_1998 = [
    PeriodicTerms([0,  0,  0,  1], [ 5.128122]),
    PeriodicTerms([0,  0,  1,  1], [ 0.280602]),
    PeriodicTerms([0,  0,  1, -1], [ 0.277693]),
    PeriodicTerms([2,  0,  0, -1], [ 0.173237]),
    PeriodicTerms([2,  0, -1,  1], [ 0.055413]),
    PeriodicTerms([2,  0, -1, -1], [ 0.046271]),
    PeriodicTerms([2,  0,  0,  1], [ 0.032573]),
    PeriodicTerms([0,  0,  2,  1], [ 0.017198]),
    PeriodicTerms([2,  0,  1, -1], [ 0.009266]),
    PeriodicTerms([0,  0,  2, -1], [ 0.008822]),
    PeriodicTerms([2, -1,  0, -1], [ 0.008216]),
    PeriodicTerms([2,  0, -2, -1], [ 0.004324]),
    PeriodicTerms([2,  0,  1,  1], [ 0.004200]),
    PeriodicTerms([2,  1,  0, -1], [-0.003359]),
    PeriodicTerms([2, -1, -1,  1], [ 0.002463]),
    PeriodicTerms([2, -1,  0,  1], [ 0.002211]),
    PeriodicTerms([2, -1, -1, -1], [ 0.002065]),
    PeriodicTerms([0,  1, -1, -1], [-0.001870]),
    PeriodicTerms([4,  0, -1, -1], [ 0.001828]),
    PeriodicTerms([0,  1,  0,  1], [-0.001794]),
    PeriodicTerms([0,  0,  0,  3], [-0.001749]),
    PeriodicTerms([0,  1, -1,  1], [-0.001565]),
    PeriodicTerms([1,  0,  0,  1], [-0.001491]),
    PeriodicTerms([0,  1,  1,  1], [-0.001475]),
    PeriodicTerms([0,  1,  1, -1], [-0.001410]),
    PeriodicTerms([0,  1,  0, -1], [-0.001344]),
    PeriodicTerms([1,  0,  0, -1], [-0.001335]),
    PeriodicTerms([0,  0,  3,  1], [ 0.001107]),
    PeriodicTerms([4,  0,  0, -1], [ 0.001021]),
    PeriodicTerms([4,  0, -1,  1], [ 0.000833]),
    PeriodicTerms([0,  0,  1, -3], [ 0.000777]),
    PeriodicTerms([4,  0, -2,  1], [ 0.000671]),
    PeriodicTerms([2,  0,  0, -3], [ 0.000607]),
    PeriodicTerms([2,  0,  2, -1], [ 0.000596]),
    PeriodicTerms([2, -1,  1, -1], [ 0.000491]),
    PeriodicTerms([2,  0, -2,  1], [-0.000451]),
    PeriodicTerms([0,  0,  3, -1], [ 0.000439]),
    PeriodicTerms([2,  0,  2,  1], [ 0.000422]),
    PeriodicTerms([2,  0, -3, -1], [ 0.000421]),
    PeriodicTerms([2,  1, -1,  1], [-0.000366]),
    PeriodicTerms([2,  1,  0,  1], [-0.000351]),
    PeriodicTerms([4,  0,  0,  1], [ 0.000331]),
    PeriodicTerms([2, -1,  1,  1], [ 0.000315]),
    PeriodicTerms([2, -2,  0, -1], [ 0.000302]),
    PeriodicTerms([0,  0,  1,  3], [-0.000283]),
    PeriodicTerms([2,  1,  1, -1], [-0.000229]),
    PeriodicTerms([1,  1,  0, -1], [ 0.000223]),
    PeriodicTerms([1,  1,  0,  1], [ 0.000223]),
    PeriodicTerms([0,  1, -2, -1], [-0.000220]),
    PeriodicTerms([2,  1, -1, -1], [-0.000220]),
    PeriodicTerms([1,  0,  1,  1], [-0.000185]),
    PeriodicTerms([2, -1, -2, -1], [ 0.000181]),
    PeriodicTerms([0,  1,  2,  1], [-0.000177]),
    PeriodicTerms([4,  0, -2, -1], [ 0.000176]),
    PeriodicTerms([4, -1, -1, -1], [ 0.000166]),
    PeriodicTerms([1,  0,  1, -1], [-0.000164]),
    PeriodicTerms([4,  0,  1, -1], [ 0.000132]),
    PeriodicTerms([1,  0, -1, -1], [-0.000119]),
    PeriodicTerms([4, -1,  0, -1], [ 0.000115]),
    PeriodicTerms([2, -2,  0,  1], [ 0.000107])]

@const_smatrix_from_series lr_1998_n lr_1998 n
@const_smatrix_from_series lr_1998_a lr_1998 a
@const_smatrix_from_series  b_1998_n  b_1998 n
@const_smatrix_from_series  b_1998_a  b_1998 a
