#   IAU 2006 Fukushima Precession-Nutation Model

# Fukashima-Williams angles
const γB_2006 = (-0.052928, 10.556378, 0.4932044, -0.00031238, -0.000002788, 0.0000000260)
const ϕB_2006 = (84381.412819, -46.811016, 0.0511268, 0.00053289, -0.000000440, -0.0000000176)
const ψB_2006 = (-0.041775, 5038.481484, 1.5584175,  -0.00018522, -0.000026452, -0.0000000148)
#   Mean obliquity of the ecliptic
const ϵB_2006 = (ϵ0_2006, -46.836769, -0.0001831, 0.00200340, -0.000000576, -0.0000000434)

#   Bias for secular variation of J2
const bj2_2006 =  0.4697e-6
#   Rate for secular variation of J2
const fj2_2006 = -2.7774e-6

#   CIO pole locator (?)
const cio_s_2006 = (94.00e-6, 3808.65e-6, -122.68e-6, -72574.11e-6, 27.98e-6, 15.62e-6)

#   Complementary equinox model constants
#   Table of coefficients of l, l', F, D, Ω, LVe, LE, pA, longitude
iau_2006_equinox_0_series = [
    #  1-10
    PeriodicTerms([ 0,  0,  0,  0,  1,  0,  0,  0], [-2640.73e-6,   0.39e-6 ]),
    PeriodicTerms([ 0,  0,  0,  0,  2,  0,  0,  0], [  -63.53e-6,   0.02e-6 ]),
    PeriodicTerms([ 0,  0,  2, -2,  3,  0,  0,  0], [  -11.75e-6,  -0.01e-6 ]),
    PeriodicTerms([ 0,  0,  2, -2,  1,  0,  0,  0], [  -11.21e-6,  -0.01e-6 ]),
    PeriodicTerms([ 0,  0,  2, -2,  2,  0,  0,  0], [    4.57e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  0,  2,  0,  3,  0,  0,  0], [   -2.02e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  0,  2,  0,  1,  0,  0,  0], [   -1.98e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  0,  0,  0,  3,  0,  0,  0], [    1.72e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  1,  0,  0,  1,  0,  0,  0], [    1.41e-6,   0.01e-6 ]),
    PeriodicTerms([ 0,  1,  0,  0, -1,  0,  0,  0], [    1.26e-6,   0.01e-6 ]),
    # 11-20
    PeriodicTerms([ 1,  0,  0,  0, -1,  0,  0,  0], [    0.63e-6,   0.00e-6 ]),
    PeriodicTerms([ 1,  0,  0,  0,  1,  0,  0,  0], [    0.63e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  1,  2, -2,  3,  0,  0,  0], [   -0.46e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  1,  2, -2,  1,  0,  0,  0], [   -0.45e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  0,  4, -4,  4,  0,  0,  0], [   -0.36e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  0,  1, -1,  1, -8, 12,  0], [    0.24e-6,   0.12e-6 ]),
    PeriodicTerms([ 0,  0,  2,  0,  0,  0,  0,  0], [   -0.32e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  0,  2,  0,  2,  0,  0,  0], [   -0.28e-6,   0.00e-6 ]),
    PeriodicTerms([ 1,  0,  2,  0,  3,  0,  0,  0], [   -0.27e-6,   0.00e-6 ]),
    PeriodicTerms([ 1,  0,  2,  0,  1,  0,  0,  0], [   -0.26e-6,   0.00e-6 ]),
    # 21-30
    PeriodicTerms([ 0,  0,  2, -2,  0,  0,  0,  0], [    0.21e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  1, -2,  2, -3,  0,  0,  0], [   -0.19e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  1, -2,  2, -1,  0,  0,  0], [   -0.18e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  0,  0,  0,  0,  8,-13, -1], [    0.10e-6,  -0.05e-6 ]),
    PeriodicTerms([ 0,  0,  0,  2,  0,  0,  0,  0], [   -0.15e-6,   0.00e-6 ]),
    PeriodicTerms([ 2,  0, -2,  0, -1,  0,  0,  0], [    0.14e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  1,  2, -2,  2,  0,  0,  0], [    0.14e-6,   0.00e-6 ]),
    PeriodicTerms([ 1,  0,  0, -2,  1,  0,  0,  0], [   -0.14e-6,   0.00e-6 ]),
    PeriodicTerms([ 1,  0,  0, -2, -1,  0,  0,  0], [   -0.14e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  0,  4, -2,  4,  0,  0,  0], [   -0.13e-6,   0.00e-6 ]),
    # 31-33
    PeriodicTerms([ 0,  0,  2, -2,  4,  0,  0,  0], [    0.11e-6,   0.00e-6 ]),
    PeriodicTerms([ 1,  0, -2,  0, -3,  0,  0,  0], [   -0.11e-6,   0.00e-6 ]),
    PeriodicTerms([ 1,  0, -2,  0, -1,  0,  0,  0], [   -0.11e-6,   0.00e-6 ])]

iau_2006_equinox_1_series = [
    #  1 - 3
    PeriodicTerms([ 0,  0,  0,  0,  2,  0,  0,  0], [   -0.07e-6,   3.57e-6 ]),
    PeriodicTerms([ 0,  0,  0,  0,  1,  0,  0,  0], [    1.73e-6,  -0.03e-6 ]),
    PeriodicTerms([ 0,  0,  2, -2,  3,  0,  0,  0], [    0.00e-6,   0.48e-6 ])]

iau_2006_equinox_2_series = [
    #  1-10
    PeriodicTerms([ 0,  0,  0,  0,  1,  0,  0,  0], [  743.52e-6,  -0.17e-6 ]),
    PeriodicTerms([ 0,  0,  2, -2,  2,  0,  0,  0], [   56.91e-6,   0.06e-6 ]),
    PeriodicTerms([ 0,  0,  2,  0,  2,  0,  0,  0], [    9.84e-6,  -0.01e-6 ]),
    PeriodicTerms([ 0,  0,  0,  0,  2,  0,  0,  0], [   -8.85e-6,   0.01e-6 ]),
    PeriodicTerms([ 0,  1,  0,  0,  0,  0,  0,  0], [   -6.38e-6,  -0.05e-6 ]),
    PeriodicTerms([ 1,  0,  0,  0,  0,  0,  0,  0], [   -3.07e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  1,  2, -2,  2,  0,  0,  0], [    2.23e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  0,  2,  0,  1,  0,  0,  0], [    1.67e-6,   0.00e-6 ]),
    PeriodicTerms([ 1,  0,  2,  0,  2,  0,  0,  0], [    1.30e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  1, -2,  2, -2,  0,  0,  0], [    0.93e-6,   0.00e-6 ]),
    # 11-20
    PeriodicTerms([ 1,  0,  0, -2,  0,  0,  0,  0], [    0.68e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  0,  2, -2,  1,  0,  0,  0], [   -0.55e-6,   0.00e-6 ]),
    PeriodicTerms([ 1,  0, -2,  0, -2,  0,  0,  0], [    0.53e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  0,  0,  2,  0,  0,  0,  0], [   -0.27e-6,   0.00e-6 ]),
    PeriodicTerms([ 1,  0,  0,  0,  1,  0,  0,  0], [   -0.27e-6,   0.00e-6 ]),
    PeriodicTerms([ 1,  0, -2, -2, -2,  0,  0,  0], [   -0.26e-6,   0.00e-6 ]),
    PeriodicTerms([ 1,  0,  0,  0, -1,  0,  0,  0], [   -0.25e-6,   0.00e-6 ]),
    PeriodicTerms([ 1,  0,  2,  0,  1,  0,  0,  0], [    0.22e-6,   0.00e-6 ]),
    PeriodicTerms([ 2,  0,  0, -2,  0,  0,  0,  0], [   -0.21e-6,   0.00e-6 ]),
    PeriodicTerms([ 2,  0, -2,  0, -1,  0,  0,  0], [    0.20e-6,   0.00e-6 ]),
    # 21-25
    PeriodicTerms([ 0,  0,  2,  2,  2,  0,  0,  0], [    0.17e-6,   0.00e-6 ]),
    PeriodicTerms([ 2,  0,  2,  0,  2,  0,  0,  0], [    0.13e-6,   0.00e-6 ]),
    PeriodicTerms([ 2,  0,  0,  0,  0,  0,  0,  0], [   -0.13e-6,   0.00e-6 ]),
    PeriodicTerms([ 1,  0,  2, -2,  2,  0,  0,  0], [   -0.12e-6,   0.00e-6 ]),
    PeriodicTerms([ 0,  0,  2,  0,  0,  0,  0,  0], [   -0.11e-6,   0.00e-6 ])]

iau_2006_equinox_3_series = [
    #  1- 4
    PeriodicTerms([ 0,  0,  0,  0,  1,  0,  0,  0], [    0.30e-6, -23.42e-6 ]),
    PeriodicTerms([ 0,  0,  2, -2,  2,  0,  0,  0], [   -0.03e-6,  -1.46e-6 ]),
    PeriodicTerms([ 0,  0,  2,  0,  2,  0,  0,  0], [   -0.01e-6,  -0.25e-6 ]),
    PeriodicTerms([ 0,  0,  0,  0,  2,  0,  0,  0], [    0.00e-6,   0.23e-6 ])]

iau_2006_equinox_4_series = [
    #  1- 1
    PeriodicTerms([ 0,  0,  0,  0,  1,  0,  0,  0], [   -0.26e-6,  -0.01e-6 ])]


@const_smatrix_from_series ϕ0_2006_equinox iau_2006_equinox_0_series n
@const_smatrix_from_series a0_2006_equinox iau_2006_equinox_0_series a
@const_smatrix_from_series ϕ1_2006_equinox iau_2006_equinox_1_series n
@const_smatrix_from_series a1_2006_equinox iau_2006_equinox_1_series a
@const_smatrix_from_series ϕ2_2006_equinox iau_2006_equinox_2_series n
@const_smatrix_from_series a2_2006_equinox iau_2006_equinox_2_series a
@const_smatrix_from_series ϕ3_2006_equinox iau_2006_equinox_3_series n
@const_smatrix_from_series a3_2006_equinox iau_2006_equinox_3_series a
@const_smatrix_from_series ϕ4_2006_equinox iau_2006_equinox_4_series n
@const_smatrix_from_series a4_2006_equinox iau_2006_equinox_4_series a
