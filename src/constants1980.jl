#   IAU 1976 model

#  sin and cos of J2000.0 mean obliquity (IAU 1976).
const COSEPS = 0.9174820620691818
const SINEPS = 0.3977771559319137

#   Precession model constants (Lieske 1979)
#   Polynomial constants for reference epoch J2000.0
const ζT_1976 = (2306.2181,  1.39656, -0.000139)
const θT_1976 = (2004.3109, -0.85330, -0.000217)
#   Polynomial constants for arbitrary epochs
#   Mean equator of epoch
const ζA_1976 = ( 0.30188, -0.000344, 0.017998)
#   Mean equator of date
const zA_1976 = ( 1.09468,  0.000066, 0.018203)
#   Difference between mean equator of epoch and date
const θA_1976 = (-0.42665, -0.000217, -0.041833)

#   IAU 1980 model

#   Initial obliquity of the ecliptic 
const ϵ0_1980  = 84381.448
#   Mean obliquity of the ecliptic (from SOFA)
const ϵ_1980   = (ϵ0_1980, -46.8150, -0.00059, 0.001813)
#   Mean longitude of the Sun minus mean longitude of Sun's perigee (from SOFA)
const l1_1980t = 99.0
const l1_1980  = (1287099.804, 1292581.224, -0.577, -0.012)
#   Mean elogation of the Moon from the Sun (from SOFA)
const D_1980t  = 1236.0
const D_1980   = (1072261.307, 1105601.328, -6.891, 0.019)
#   Mean longitude of the Moon minus mean longitude of Moon's perigee (from SOFA)
const l0_1980t = 1325.0
const l0_1980  = (485866.733, 715922.633, 31.310, 0.064)
#   Mean longitude of the Moon minus mean longitude of Moon's node (from SOFA)
const F_1980t  = 1342.0
const F_1980   = (335778.877, 295263.137, -13.257, 0.011)
#   Mean Longitude of the Moon ascending node on the ecliptic,
#   measured from the mean equinox of date (from SOFA)
const Ω_1980t  = -5.0
const Ω_1980   = (450160.280, -482890.539, 7.455, 0.008)

#   Nutation series model constants
#   (see Explanatory Supplement, Table 3.222.1)
#   Table of coefficents of 1, l', F, D, Ω, longitude, obliquity
const iau_1980_nutation_series = [
    #  1-10
    PeriodicTerms([ 0,  0,  0,  0,  1], [ -171996.0, -174.2,  92025.0,    8.9]),
    PeriodicTerms([ 0,  0,  0,  0,  2], [    2062.0,    0.2,   -895.0,    0.5]),
    PeriodicTerms([-2,  0,  2,  0,  1], [      46.0,    0.0,    -24.0,    0.0]),
    PeriodicTerms([ 2,  0, -2,  0,  0], [      11.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([-2,  0,  2,  0,  2], [      -3.0,    0.0,      1.0,    0.0]),
    PeriodicTerms([ 1, -1,  0, -1,  0], [      -3.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0, -2,  2, -2,  1], [      -2.0,    0.0,      1.0,    0.0]),
    PeriodicTerms([ 2,  0, -2,  0,  1], [       1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0,  0,  2, -2,  2], [  -13187.0,   -1.6,   5736.0,   -3.1]),
    PeriodicTerms([ 0,  1,  0,  0,  0], [    1426.0,   -3.4,     54.0,   -0.1]),
    # 11-20
    PeriodicTerms([ 0,  1,  2, -2,  2], [    -517.0,    1.2,    224.0,   -0.6]),
    PeriodicTerms([ 0, -1,  2, -2,  2], [     217.0,   -0.5,    -95.0,    0.3]),
    PeriodicTerms([ 0,  0,  2, -2,  1], [     129.0,    0.1,    -70.0,    0.0]),
    PeriodicTerms([ 2,  0,  0, -2,  0], [      48.0,    0.0,      1.0,    0.0]),
    PeriodicTerms([ 0,  0,  2, -2,  0], [     -22.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0,  2,  0,  0,  0], [      17.0,   -0.1,      0.0,    0.0]),
    PeriodicTerms([ 0,  1,  0,  0,  1], [     -15.0,    0.0,      9.0,    0.0]),
    PeriodicTerms([ 0,  2,  2, -2,  2], [     -16.0,    0.1,      7.0,    0.0]),
    PeriodicTerms([ 0, -1,  0,  0,  1], [     -12.0,    0.0,      6.0,    0.0]),
    PeriodicTerms([-2,  0,  0,  2,  1], [      -6.0,    0.0,      3.0,    0.0]),
    # 21-30
    PeriodicTerms([ 0, -1,  2, -2,  1], [      -5.0,    0.0,      3.0,    0.0]),
    PeriodicTerms([ 2,  0,  0, -2,  1], [       4.0,    0.0,     -2.0,    0.0]),
    PeriodicTerms([ 0,  1,  2, -2,  1], [       4.0,    0.0,     -2.0,    0.0]),
    PeriodicTerms([ 1,  0,  0, -1,  0], [      -4.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 2,  1,  0, -2,  0], [       1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0,  0, -2,  2,  1], [       1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0,  1, -2,  2,  0], [      -1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0,  1,  0,  0,  2], [       1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([-1,  0,  0,  1,  1], [       1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0,  1,  2, -2,  0], [      -1.0,    0.0,      0.0,    0.0]),
    # 31-40
    PeriodicTerms([ 0,  0,  2,  0,  2], [   -2274.0,   -0.2,    977.0,   -0.5]),
    PeriodicTerms([ 1,  0,  0,  0,  0], [     712.0,    0.1,     -7.0,    0.0]),
    PeriodicTerms([ 0,  0,  2,  0,  1], [    -386.0,   -0.4,    200.0,    0.0]),
    PeriodicTerms([ 1,  0,  2,  0,  2], [    -301.0,    0.0,    129.0,   -0.1]),
    PeriodicTerms([ 1,  0,  0, -2,  0], [    -158.0,    0.0,     -1.0,    0.0]),
    PeriodicTerms([-1,  0,  2,  0,  2], [     123.0,    0.0,    -53.0,    0.0]),
    PeriodicTerms([ 0,  0,  0,  2,  0], [      63.0,    0.0,     -2.0,    0.0]),
    PeriodicTerms([ 1,  0,  0,  0,  1], [      63.0,    0.1,    -33.0,    0.0]),
    PeriodicTerms([-1,  0,  0,  0,  1], [     -58.0,   -0.1,     32.0,    0.0]),
    PeriodicTerms([-1,  0,  2,  2,  2], [     -59.0,    0.0,     26.0,    0.0]),
    # 41-50
    PeriodicTerms([ 1,  0,  2,  0,  1], [     -51.0,    0.0,     27.0,    0.0]),
    PeriodicTerms([ 0,  0,  2,  2,  2], [     -38.0,    0.0,     16.0,    0.0]),
    PeriodicTerms([ 2,  0,  0,  0,  0], [      29.0,    0.0,     -1.0,    0.0]),
    PeriodicTerms([ 1,  0,  2, -2,  2], [      29.0,    0.0,    -12.0,    0.0]),
    PeriodicTerms([ 2,  0,  2,  0,  2], [     -31.0,    0.0,     13.0,    0.0]),
    PeriodicTerms([ 0,  0,  2,  0,  0], [      26.0,    0.0,     -1.0,    0.0]),
    PeriodicTerms([-1,  0,  2,  0,  1], [      21.0,    0.0,    -10.0,    0.0]),
    PeriodicTerms([-1,  0,  0,  2,  1], [      16.0,    0.0,     -8.0,    0.0]),
    PeriodicTerms([ 1,  0,  0, -2,  1], [     -13.0,    0.0,      7.0,    0.0]),
    PeriodicTerms([-1,  0,  2,  2,  1], [     -10.0,    0.0,      5.0,    0.0]),
    # 51-60
    PeriodicTerms([ 1,  1,  0, -2,  0], [      -7.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0,  1,  2,  0,  2], [       7.0,    0.0,     -3.0,    0.0]),
    PeriodicTerms([ 0, -1,  2,  0,  2], [      -7.0,    0.0,      3.0,    0.0]),
    PeriodicTerms([ 1,  0,  2,  2,  2], [      -8.0,    0.0,      3.0,    0.0]),
    PeriodicTerms([ 1,  0,  0,  2,  0], [       6.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 2,  0,  2, -2,  2], [       6.0,    0.0,     -3.0,    0.0]),
    PeriodicTerms([ 0,  0,  0,  2,  1], [      -6.0,    0.0,      3.0,    0.0]),
    PeriodicTerms([ 0,  0,  2,  2,  1], [      -7.0,    0.0,      3.0,    0.0]),
    PeriodicTerms([ 1,  0,  2, -2,  1], [       6.0,    0.0,     -3.0,    0.0]),
    PeriodicTerms([ 0,  0,  0, -2,  1], [      -5.0,    0.0,      3.0,    0.0]),
    # 61-70
    PeriodicTerms([ 1, -1,  0,  0,  0], [       5.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 2,  0,  2,  0,  1], [      -5.0,    0.0,      3.0,    0.0]),
    PeriodicTerms([ 0,  1,  0, -2,  0], [      -4.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 1,  0, -2,  0,  0], [       4.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0,  0,  0,  1,  0], [      -4.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 1,  1,  0,  0,  0], [      -3.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 1,  0,  2,  0,  0], [       3.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 1, -1,  2,  0,  2], [      -3.0,    0.0,      1.0,    0.0]),
    PeriodicTerms([-1, -1,  2,  2,  2], [      -3.0,    0.0,      1.0,    0.0]),
    PeriodicTerms([-2,  0,  0,  0,  1], [      -2.0,    0.0,      1.0,    0.0]),
    # 71-80
    PeriodicTerms([ 3,  0,  2,  0,  2], [      -3.0,    0.0,      1.0,    0.0]),
    PeriodicTerms([ 0, -1,  2,  2,  2], [      -3.0,    0.0,      1.0,    0.0]),
    PeriodicTerms([ 1,  1,  2,  0,  2], [       2.0,    0.0,     -1.0,    0.0]),
    PeriodicTerms([-1,  0,  2, -2,  1], [      -2.0,    0.0,      1.0,    0.0]),
    PeriodicTerms([ 2,  0,  0,  0,  1], [       2.0,    0.0,     -1.0,    0.0]),
    PeriodicTerms([ 1,  0,  0,  0,  2], [      -2.0,    0.0,      1.0,    0.0]),
    PeriodicTerms([ 3,  0,  0,  0,  0], [       2.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0,  0,  2,  1,  2], [       2.0,    0.0,     -1.0,    0.0]),
    PeriodicTerms([-1,  0,  0,  0,  2], [       1.0,    0.0,     -1.0,    0.0]),
    PeriodicTerms([ 1,  0,  0, -4,  0], [      -1.0,    0.0,      0.0,    0.0]),
    # 81-90
    PeriodicTerms([-2,  0,  2,  2,  2], [       1.0,    0.0,     -1.0,    0.0]),
    PeriodicTerms([-1,  0,  2,  4,  2], [      -2.0,    0.0,      1.0,    0.0]),
    PeriodicTerms([ 2,  0,  0, -4,  0], [      -1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 1,  1,  2, -2,  2], [       1.0,    0.0,     -1.0,    0.0]),
    PeriodicTerms([ 1,  0,  2,  2,  1], [      -1.0,    0.0,      1.0,    0.0]),
    PeriodicTerms([-2,  0,  2,  4,  2], [      -1.0,    0.0,      1.0,    0.0]),
    PeriodicTerms([-1,  0,  4,  0,  2], [       1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 1, -1,  0, -2,  0], [       1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 2,  0,  2, -2,  1], [       1.0,    0.0,     -1.0,    0.0]),
    PeriodicTerms([ 2,  0,  2,  2,  2], [      -1.0,    0.0,      0.0,    0.0]),
    # 91-100
    PeriodicTerms([ 1,  0,  0,  2,  1], [      -1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0,  0,  4, -2,  2], [       1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 3,  0,  2, -2,  2], [       1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 1,  0,  2, -2,  0], [      -1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0,  1,  2,  0,  1], [       1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([-1, -1,  0,  2,  1], [       1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0,  0, -2,  0,  1], [      -1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0,  0,  2, -1,  2], [      -1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0,  1,  0,  2,  0], [      -1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 1,  0, -2, -2,  0], [      -1.0,    0.0,      0.0,    0.0]),
    # 101-106
    PeriodicTerms([ 0, -1,  2,  0,  1], [      -1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 1,  1,  0, -2,  1], [      -1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 1,  0, -2,  2,  0], [      -1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 2,  0,  0,  2,  0], [       1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0,  0,  2,  4,  2], [      -1.0,    0.0,      0.0,    0.0]),
    PeriodicTerms([ 0,  1,  0,  1,  0], [       1.0,    0.0,      0.0,    0.0])]

@const_smatrix_from_series n_1980_nutation iau_1980_nutation_series n
@const_smatrix_from_series a_1980_nutation iau_1980_nutation_series a
