####    Astronomy / Astrometry    ####

"""
	ab(pnat::AbstractVector{<:AbstractFloat}, v::AbstractVector{<:AbstractFloat}, s::AbstractFloat, bm1::AbstractFloat)

Apply aberration to transform natural direction into proper direction.

# Input

 - `pnat` --  natural direction to source (unit vector).
 - `v`    -- observer barycentric velocity in units of light speed.
 - `s`    -- distance between the Sun and the observer in units of AU.
 - `bm1`  -- inverse Lorenz factor, i.e., (sqrt(1-v^2)).

# Output

 - `p`    -- proper direction to source (unit vector).

# Note

1) The algorithm is based on Expr. (7.40) in the Explanatory
   Supplement (Urban & Seidelmann 2013), but with the following
   changes:

   o  Rigorous rather than approximate normalization is applied.

   o  The gravitational potential term from Expr. (7) in Klioner
	  (2003) is added, taking into account only the Sun's
	  contribution.  This has a maximum effect of about 0.4
	  microarcsecond.

2) In almost all cases, the maximum accuracy will be limited by the
   supplied velocity.  For example, if the ERFA eraEpv00 function is
   used, errors of up to 5 microarcseconds could occur.

# References

Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to the
Astronomical Almanac, 3rd ed., University Science Books (2013).

Klioner, Sergei A., "A practical relativistic model for micro-
arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).
"""
function ab(pnat::V, v::W, s::F, bm1::F) where
{V <: AbstractVector{<:AbstractFloat}, W <: AbstractVector{<:AbstractFloat},
	F <: AbstractFloat}
	p = bm1 .* pnat .+ (1.0 + sum(pnat .* v)/(1.0 + bm1)) .* v .+
		SCHWARZRADIUS/s .* (v .- sum(pnat .* v) .* pnat)
	p ./ norm(p)
end

"""
	apcg(day1::AbstractFloat, day2::AbstractFloat, ebpv::AbstractVector{<:AbstractVector{<:AbstractFloat}},
		 ehp::AbstractVector{<:AbstractFloat})

For a geocentric observer, prepare star-independent astrometry
paramenters for transformations between ICRS and GCRS coordinates. The
Earth ephemeris is supplied by the caller.

The parameters produced by this function are required in the parallax,
light deflection and aberration parts of the astrometric
transformation chain.
 
# Input

 - `day1`   -- TDB as a 2-part...
 - `day2`   -- ...Julian Date (Note 1)
 - `ebpv`   -- Earth barycentric pos/vel (AU, AU/day)
 - `ehp`    -- Earth heliocentric position (AU)

# Output

 - `astrom` -- star-independent astrometry parameters:

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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

2) All the vectors are with respect to BCRS axes.

3) This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

		functions         observer        transformation

	   apcg apcg13    geocentric      ICRS <-> GCRS
	   apci apci13    terrestrial     ICRS <-> CIRS
	   apco apco13    terrestrial     ICRS <-> observed
	   apcs apcs13    space           ICRS <-> GCRS
	   aper aper13    terrestrial     update Earth rotation
	   apio apio13    terrestrial     CIRS <-> observed

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion, parallax,
   light deflection, and aberration.  From GCRS to CIRS comprises
   frame bias and precession-nutation.  From CIRS to observed takes
   account of Earth rotation, polar motion, diurnal aberration and
   parallax (unless subsumed into the ICRS <-> GCRS transformation),
   and atmospheric refraction.

4) The context structure astrom produced by this function is used by
   eraAtciq* and eraAticq*.
"""
function apcg(day1::F, day2::F, ebpv::W, ehp::V) where
{F <: AbstractFloat, V <: AbstractVector{<:AbstractFloat},
	W <: AbstractVecOrMat{<:AbstractVector{<:AbstractFloat}}}

	#  Compute the star-independent astrometry parameters.
	apcs(day1, day2, SVector(SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0)), ebpv, ehp)
end

"""
	apcg13(day1::AbstractFloat, day2::AbstractFloat)

For a geocentric observer, prepare star-independent astrometry
parameters for transformations between ICRS and GCRS coordinates.  The
caller supplies the date, and ERFA models are used to predict the
Earth ephemeris.

The parameters produced by this function are required in the parallax,
light deflection and aberration parts of the astrometric
transformation chain.

# Input

 - `day1`   -- TDB as a 2-part...
 - `day2`   -- ...Julian Date (Note 1)

# Output

 - `astrom` -- star-independent astrometry parameters:

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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

2) All the vectors are with respect to BCRS axes.

3) In cases where the caller wishes to supply his own Earth
   ephemeris, the function eraApcg can be used instead of the present
   function.

4) This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

	  functions       observer        transformation

	 apcg apcg13    geocentric      ICRS <-> GCRS
	 apci apci13    terrestrial     ICRS <-> CIRS
	 apco apco13    terrestrial     ICRS <-> observed
	 apcs apcs13    space           ICRS <-> GCRS
	 aper aper13    terrestrial     update Earth rotation
	 apio apio13    terrestrial     CIRS <-> observed

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion, parallax,
   light deflection, and aberration.  From GCRS to CIRS comprises
   frame bias and precession-nutation.  From CIRS to observed takes
   account of Earth rotation, polar motion, diurnal aberration and
   parallax (unless subsumed into the ICRS <-> GCRS transformation),
   and atmospheric refraction.

5) The context structure astrom produced by this function is used by
   atciq* and aticq*.
"""
function apcg13(day1::AbstractFloat, day2::AbstractFloat)
	#  Earth barycentric and heliocentric position and velocity (AU, AU/day).
	ehpv, ebpv = epv00(day1, day2)
	#  Compute the star-independent astrometry parameters.
	apcg(day1, day2, ebpv, ehpv[1])
end

"""
	apci(day1::AbstractFloat, day2::AbstractFloat, ebpv::AbstractVector{<:AbstractVector{<:AbstractFloat}},
		 ehp::AbstractVector{<:AbstractFloat}, x::AbstractFloat, y::AbstractFloat, s::AbstractFloat)

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between ICRS and geocentric CIRS
coordinates.  The Earth ephemeris and CIP/CIO are supplied by the
caller.

The parameters produced by this function are required in the parallax,
light deflection, aberration, and bias-precession-nutation parts of
the astrometric transformation chain.

# Input

 - `day1`   -- TDB as a 2-part...
 - `day2`   -- ...Julian Date (Note 1)
 - `ebpv`   -- Earth barycentric position/velocity (AU, AU/day)
 - `ehp`    -- Earth heliocentric position (AU)
 - `x, y`   -- CIP X,Y (components of unit vector)
 - `s`      -- the CIO locator s (radians)

# Output

 - `astrom` -- star-independent astrometry parameters

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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

2) All the vectors are with respect to BCRS axes.

3) In cases where the caller does not wish to provide the Earth
   ephemeris and CIP/CIO, the function eraApci13 can be used instead
   of the present function.  This computes the required quantities
   using other ERFA functions.

4) This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

	  functions       observer        transformation

	 apcg apcg13    geocentric      ICRS <-> GCRS
	 apci apci13    terrestrial     ICRS <-> CIRS
	 apco apco13    terrestrial     ICRS <-> observed
	 apcs apcs13    space           ICRS <-> GCRS
	 aper aper13    terrestrial     update Earth rotation
	 apio apio13    terrestrial     CIRS <-> observed

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion, parallax,
   light deflection, and aberration.  From GCRS to CIRS comprises
   frame bias and precession-nutation.  From CIRS to observed takes
   account of Earth rotation, polar motion, diurnal aberration and
   parallax (unless subsumed into the ICRS <-> GCRS transformation),
   and atmospheric refraction.

5) The context structure astrom produced by this function is used by
   atciq and aticq.
"""
function apci(day1::F, day2::F, ebpv::W, ehp::V, x::F, y::F, s::F) where
{F <: AbstractFloat, V <: AbstractVector{<:AbstractFloat},
	W <: AbstractVecOrMat{<:AbstractVector{<:AbstractFloat}}}

	#  Star-independent astrometry parameters for geocenter and CIO based
	#  bias-precession-nutation matrix.
	p = apcg(day1, day2, ebpv, ehp)
	Astrom(p.pmt, p.eb, p.eh, p.em, p.v, p.bm1, c2ixys(x, y, s))
end

"""
	apci13(day1::AbstractFloat, day2::AbstractFloat)

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between ICRS and geocentric CIRS
coordinates.  The caller supplies the date, and ERFA models are used
to predict the Earth ephemeris and CIP/CIO.

The parameters produced by this function are required in the parallax,
light deflection, aberration, and bias-precession-nutation parts of
the astrometric transformation chain.

# Input

 - `date1`  -- TDB as a 2-part...
 - `date2`  -- ...Julian Date (Note 1)

# Output

 - `astrom` -- star-independent astrometry parameters:

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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

2) All the vectors are with respect to BCRS axes.

3) In cases where the caller wishes to supply his own Earth
   ephemeris and CIP/CIO, the function eraApci can be used instead
   of the present function.

4) This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

	  functions       observer        transformation

	 apcg apcg13    geocentric      ICRS <-> GCRS
	 apci apci13    terrestrial     ICRS <-> CIRS
	 apco apco13    terrestrial     ICRS <-> observed
	 apcs apcs13    space           ICRS <-> GCRS
	 aper aper13    terrestrial     update Earth rotation
	 apio apio13    terrestrial     CIRS <-> observed

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion, parallax,
   light deflection, and aberration.  From GCRS to CIRS comprises
   frame bias and precession-nutation.  From CIRS to observed takes
   account of Earth rotation, polar motion, diurnal aberration and
   parallax (unless subsumed into the ICRS <-> GCRS transformation),
   and atmospheric refraction.

5) The context structure astrom produced by this function is used by
   atciq and aticq.
"""
function apci13(day1::AbstractFloat, day2::AbstractFloat)
	#  Earth barycentric and heliocentric position and velocity (AU, AU/day).
	ehpv, ebpv = epv00(day1, day2)
	#  Form the equinox based bpn matrix, IAU 2006/2000A
	bpn = pnm06a(day1, day2)
	#  Extract the CIP x and y
	x, y = bpn2xy(bpn)
	#  Obtain the CIO locator
	s = s06(day1, day2, x, y)
	#  Compute the star-independent astrometry parameters.
	astrom = apci(day1, day2, ebpv, ehpv[1], x, y, s)
	#  Return the star-independent parameters and equation of origin.
	(astrom, eors(bpn, s))
end

"""
	apco(day1::AbstractFloat, day2::AbstractFloat, ebpv::AbstractVector{<:AbstractVector{<:AbstractFloat}},
		 ehp::AbstractVector{<:AbstractFloat}, x::AbstractFloat, y::AbstractFloat, s::AbstractFloat,
		 θ::AbstractFloat, elong::AbstractFloat, ϕ::AbstractFloat, hm::AbstractFloat,
		 xp::AbstractFloat, yp::AbstractFloat, sp::AbstractFloat, refa::AbstractFloat, refb::AbstractFloat)

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between ICRS and observed coordinates.
The caller supplies the Earth ephemeris, the Earth rotation
information and the refraction constants as well as the site
coordinates.

# Input

 - `date1`  -- TDB as a 2-part...
 - `date2`  -- ...Julian Date (Note 1)
 - `ebpv`   -- Earth barycentric PV (AU, AU/day, Note 2)
 - `ehp`    -- Earth heliocentric P (AU, Note 2)
 - `x, y`   -- CIP X,Y (components of unit vector)
 - `s`      -- the CIO locator s (radians)
 - `theta`  -- Earth rotation angle (radians)
 - `elong`  -- longitude (radians, east +ve, Note 3)
 - `phi`    -- latitude (geodetic, radians, Note 3)
 - `hm`     -- height above ellipsoid (m, geodetic, Note 3)
 - `xp, yp` -- polar motion coordinates (radians, Note 4)
 - `sp`     -- the TIO locator s' (radians, Note 4)
 - `refa`   -- refraction constant A (radians, Note 5)
 - `refb`   -- refraction constant B (radians, Note 5)

# Output

 - `astrom` -- star-independent astrometry parameters:

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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

2) The vectors eb, eh, and all the astrom vectors, are with respect
   to BCRS axes.

3) The geographical coordinates are with respect to the ERFA_WGS84
   reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN CONVENTION:
   the longitude required by the present function is right-handed,
   i.e. east-positive, in accordance with geographical convention.

   The adjusted longitude stored in the astrom array takes into
   account the TIO locator and polar motion.

4) xp and yp are the coordinates (in radians) of the Celestial
   Intermediate Pole with respect to the International Terrestrial
   Reference System (see IERS Conventions), measured along the
   meridians 0 and 90 deg west respectively.  sp is the TIO locator
   s', in radians, which positions the Terrestrial Intermediate Origin
   on the equator.  For many applications, xp, yp and (especially) sp
   can be set to zero.

   Internally, the polar motion is stored in a form rotated onto the
   local meridian.

5) The refraction constants refa and refb are for use in a dZ =
   A*tan(Z)+B*tan^3(Z) model, where Z is the observed (i.e. refracted)
   zenith distance and dZ is the amount of refraction.

6) It is advisable to take great care with units, as even unlikely
   values of the input parameters are accepted and processed in
   accordance with the models used.

7) In cases where the caller does not wish to provide the Earth
   Ephemeris, the Earth rotation information and refraction constants,
   the function eraApco13 can be used instead of the present function.
   This starts from UTC and weather readings etc.  and computes
   suitable values using other ERFA functions.

8) This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

		functions         observer        transformation

	 eraApcg eraApcg13    geocentric      ICRS <-> GCRS
	 eraApci eraApci13    terrestrial     ICRS <-> CIRS
	 eraApco eraApco13    terrestrial     ICRS <-> observed
	 eraApcs eraApcs13    space           ICRS <-> GCRS
	 eraAper eraAper13    terrestrial     update Earth rotation
	 eraApio eraApio13    terrestrial     CIRS <-> observed

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion, parallax,
   light deflection, and aberration.  From GCRS to CIRS comprises
   frame bias and precession-nutation.  From CIRS to observed takes
   account of Earth rotation, polar motion, diurnal aberration and
   parallax (unless subsumed into the ICRS <-> GCRS transformation),
   and atmospheric refraction.

9) The context structure astrom produced by this function is used by
   eraAtioq, eraAtoiq, eraAtciq* and eraAticq*.
"""
function apco(day1::F, day2::F, ebpv::W, ehp::V, x::F, y::F, s::F, θ::F,
	elong::F, ϕ::F, hm::F, xp::F, yp::F, sp::F, refa::F, refb::F) where
	{F <: AbstractFloat, V <: AbstractVector{<:AbstractFloat},
	W <: AbstractVecOrMat{<:AbstractVector{<:AbstractFloat}}}

	#  Form the rotation matrix, CIRS to apparent (HA, Dec).
	r = Rz(elong) * Rx(-yp) * Ry(-xp) * Rz(θ+sp)
	#  Solve for the local Earth rotation angle.
	eral = r[1, 1] != 0.0 || r[1, 2] != 0.0 ? atan(r[1, 2], r[1, 1]) : 0.0
	#  Solve for the polar motion (x, y) with respect to local meridian.
	xpl = atan(r[1, 3], norm(r[1, 1:2]))
	ypl = r[2, 3] != 0.0 || r[3, 3] != 0.0 ? -atan(r[2, 3], r[3, 3]) : 0.0
	#  Adjusted longitude.
	along = rem2pi(eral - θ, RoundNearest)
	#  Functions of latitude.
	sphi, cphi = sincos(ϕ)
	#  CIO based bpn matrix
	r = c2ixys(x, y, s)
	#  Observer's geocentric position and velocity (m, m/s, CIRS) and
	#  rotate into GCRS.
	pv = trxpv(r, pvtob(elong, ϕ, hm, xp, yp, sp, θ))
	#  ICRS <-> GCRS parameters.
	a = apcs(day1, day2, pv, ebpv, ehp)
	Astrom(a.pmt, a.eb, a.eh, a.em, a.v, a.bm1, r, along, a.phi, xpl, ypl,
		sphi, cphi, 0.0, eral, refa, refb)
end

"""
	apco13(day1::AbstractFloat, day2::AbstractFloat, dut1::AbstractFloat, elong::AbstractFloat,
		   phi::AbstractFloat, hm::AbstractFloat, xp::AbstractFloat, yp::AbstractFloat,
		   phpa::AbstractFloat, tc::AbstractFloat, rh::AbstractFloat, wl::AbstractFloat)

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between ICRS and observed coordinates.
The caller supplies UTC, site coordinates, ambient air conditions and
observing wavelength, and ERFA models are used to obtain the Earth
ephemeris, CIP/CIO and refraction constants.

The parameters produced by this function are required in the parallax,
light deflection, aberration, and bias-precession-nutation parts of
the ICRS/CIRS transformations.

# Input

 - `utc1`   -- UTC as a 2-part...
 - `utc2`   -- ...quasi Julian Date (Notes 1,2)
 - `dut1`   -- UT1-UTC (seconds, Note 3)
 - `elong`  -- longitude (radians, east +ve, Note 4)
 - `phi`    -- latitude (geodetic, radians, Note 4)
 - `hm`     -- height above ellipsoid (m, geodetic, Notes 4,6)
 - `xp, yp` -- polar motion coordinates (radians, Note 5)
 - `phpa`   -- pressure at the observer (hPa = mB, Note 6)
 - `tc`     -- ambient temperature at the observer (deg C)
 - `rh`     -- relative humidity at the observer (range 0-1)
 - `wl`     -- wavelength (micrometers, Note 7)

# Output

 - astrom`  -- star-independent astrometry parameters:

# Note

1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
   convenient way between the two arguments, for example where utc1 is
   the Julian Day Number and utc2 is the fraction of a day.

   However, JD cannot unambiguously represent UTC during a leap second
   unless special measures are taken.  The convention in the present
   function is that the JD day represents UTC days whether the length
   is 86399, 86400 or 86401 SI seconds.

   Applications should use the function eraDtf2d to convert from
   calendar date and time of day into 2-part quasi Julian Date, as it
   implements the leap-second-ambiguity convention just described.

2) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale or that are too far in the future to
   be trusted.  See eraDat for further details.

3) UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
   one second at the end of each positive UTC leap second, introduced
   in order to keep UT1-UTC within +/- 0.9s.  n.b. This practice is
   under review, and in the future UT1-UTC may grow essentially
   without limit.

4) The geographical coordinates are with respect to the ERFA_WGS84
   reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN: the
   longitude required by the present function is east-positive
   (i.e. right-handed), in accordance with geographical convention.

5) The polar motion xp,yp can be obtained from IERS bulletins.  The
   values are the coordinates (in radians) of the Celestial
   Intermediate Pole with respect to the International Terrestrial
   Reference System (see IERS Conventions 2003), measured along the
   meridians 0 and 90 deg west respectively.  For many applications,
   xp and yp can be set to zero.

   Internally, the polar motion is stored in a form rotated onto the
   local meridian.

6) If hm, the height above the ellipsoid of the observing station in
   meters, is not known but phpa, the pressure in hPa (=mB), is
   available, an adequate estimate of hm can be obtained from the
   expression

		  hm = -29.3 * tsl * log ( phpa / 1013.25 );

   where tsl is the approximate sea-level air temperature in K (See
   Astrophysical Quantities, C.W.Allen, 3rd edition, section 52).
   Similarly, if the pressure phpa is not known, it can be estimated
   from the height of the observing station, hm, as follows:

		  phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

   Note, however, that the refraction is nearly proportional to the
   pressure and that an accurate phpa value is important for precise
   work.

7) The argument wl specifies the observing wavelength in micrometers.
   The transition from optical to radio is assumed to occur at 100
   micrometers (about 3000 GHz).

8) It is advisable to take great care with units, as even unlikely
   values of the input parameters are accepted and processed in
   accordance with the models used.

9) In cases where the caller wishes to supply his own Earth ephemeris,
   Earth rotation information and refraction constants, the function
   eraApco can be used instead of the present function.

10) This is one of several functions that inserts into the astrom
	structure star-independent parameters needed for the chain of
	astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

	The various functions support different classes of observer and
	portions of the transformation chain:

		functions         observer        transformation

	 eraApcg eraApcg13    geocentric      ICRS <-> GCRS
	 eraApci eraApci13    terrestrial     ICRS <-> CIRS
	 eraApco eraApco13    terrestrial     ICRS <-> observed
	 eraApcs eraApcs13    space           ICRS <-> GCRS
	 eraAper eraAper13    terrestrial     update Earth rotation
	 eraApio eraApio13    terrestrial     CIRS <-> observed

	Those with names ending in "13" use contemporary ERFA models to
	compute the various ephemerides.  The others accept ephemerides
	supplied by the caller.

	The transformation from ICRS to GCRS covers space motion,
	parallax, light deflection, and aberration.  From GCRS to CIRS
	comprises frame bias and precession-nutation.  From CIRS to
	observed takes account of Earth rotation, polar motion, diurnal
	aberration and parallax (unless subsumed into the ICRS <-> GCRS
	transformation), and atmospheric refraction.

11) The context structure astrom produced by this function is used by
	eraAtioq, eraAtoiq, eraAtciq* and eraAticq*.
"""
function apco13(day1::AbstractFloat, day2::AbstractFloat, dut1::AbstractFloat,
	elong::AbstractFloat, phi::AbstractFloat, hm::AbstractFloat, xp::AbstractFloat,
	yp::AbstractFloat, phpa::AbstractFloat, tc::AbstractFloat, rh::AbstractFloat,
	wl::AbstractFloat)
	#  UTC to other time scales.
	tt1, tt2 = taitt(utctai(day1, day2)...)
	#  Earth barycentric and heliocentric position and velocity (AU, AU/day).
	ehpv, ebpv = epv00(tt1, tt2)
	#  Form the equinox based BPN matrix (IAU 2006/2000A).
	r = pnm06a(tt1, tt2)
	#  Extract CIP x, y
	x, y = bpn2xy(r)
	#  Obtain CIO locator s and the TIO locator s'.
	s, sp = s06(tt1, tt2, x, y), sp00(tt1, tt2)
	#  Earth rotation angle
	θ = era00(utcut1(day1, day2, dut1)...)
	#  Refraction constants A and B.
	refa, refb = refco(phpa, tc, rh, wl)
	#  Compute the star-independent astrometry paramaters.
	astrom = apco(tt1, tt2, ebpv, ehpv[1], x, y, s, θ, elong, phi, hm,
		xp, yp, sp, refa, refb)
	(astrom, eors(r, s))
end

"""
	apcs(day1::AbstractFloat, day2::AbstractFloat, pv::AbstractVector{<:AbstractVector{<:AbstractFloat}},
		 ebpv::AbstractVector{<:AbstractVector{<:AbstractFloat}}, ehp::AbstractVector{<:AbstractFloat})

For an observer whose geocentric position and velocity are known,
prepare star-independent astrometry parameters for transformations
between ICRS and GCRS.  The Earth ephemeris is supplied by the caller.

The parameters produced by this function are required in the space
motion, parallax, light deflection and aberration parts of the
astrometric transformation chain.

# Input

 - `day1`   -- TDB as a 2-part...
 - `day2`   -- ...Julian Date (Note 1)
 - `pv`     -- observer's geocentric pos/vel (m, m/s)
 - `ebpv`   -- Earth barycentric PV (AU, AU/day)
 - `ehp`    -- Earth heliocentric P (AU)

# Output

 - `astrom` -- star-independent astrometry parameters:

# Note

1) The TDB date date1+date2 is a Julian Date, apportioned in any

   convenient way between the two arguments.  For example,
   JD(TDB)=2450123.7 could be expressed in any of these ways, among
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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

2) All the vectors are with respect to BCRS axes.

3) Providing separate arguments for (i) the observer's geocentric
   position and velocity and (ii) the Earth ephemeris is done for
   convenience in the geocentric, terrestrial and Earth orbit cases.
   For deep space applications it maybe more convenient to specify
   zero geocentric position and velocity and to supply the observer's
   position and velocity information directly instead of with respect
   to the Earth.  However, note the different units: m and m/s for the
   geocentric vectors, AU and AU/day for the heliocentric and
   barycentric vectors.

4) In cases where the caller does not wish to provide the Earth
   ephemeris, the function eraApcs13 can be used instead of the
   present function.  This computes the Earth ephemeris using the ERFA
   function eraEpv00.

5) This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

		functions        observer        transformation

	 Apcg eraApcg13    geocentric      ICRS <-> GCRS
	 Apci eraApci13    terrestrial     ICRS <-> CIRS
	 Apco eraApco13    terrestrial     ICRS <-> observed
	 Apcs eraApcs13    space           ICRS <-> GCRS
	 Aper eraAper13    terrestrial     update Earth rotation
	 Apio eraApio13    terrestrial     CIRS <-> observed

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion, parallax,
   light deflection, and aberration.  From GCRS to CIRS comprises
   frame bias and precession-nutation.  From CIRS to observed takes
   account of Earth rotation, polar motion, diurnal aberration and
   parallax (unless subsumed into the ICRS <-> GCRS transformation),
   and atmospheric refraction.

6) The context structure astrom produced by this function is used by
   Atciq and Aticq.
"""
function apcs(day1::F, day2::F, pv::W, ebpv::X, ehp::V) where
{F <: AbstractFloat, V <: AbstractVector{<:AbstractFloat},
	W <: AbstractVector{<:AbstractVector{<:AbstractFloat}},
	X <: AbstractVector{<:AbstractVector{<:AbstractFloat}}}
	# Time since reference epoch, years (for proper motion calculation).
	pmt = ((day1 - JD2000) + day2)/DAYPERYEAR
	# Barycentric position of observer (AU).
	eb = ebpv[1] .+ pv[1] ./ ASTRUNIT
	# Barycentric velocity (in speed of light).
	ev = (ebpv[2] .+ pv[2] ./ (ASTRUNIT/SECPERDAY)) .* (ASTRUNIT/LIGHTSPEED/SECPERDAY)
	# Heliocentric direction and distance (unit vector & AU).
	em, eh = pn(ehp .+ pv[1] ./ ASTRUNIT)
	#  Recprocal of Lorenz factor
	bm1 = sqrt(1.0 - sum(ev .* ev))
	Astrom(pmt, eb, eh, em, ev, bm1,
		SMatrix{3, 3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
end

"""
	apcs13(day1::AbstractFloat, day2::AbstractFloat, pv::AbstractVector{<:AbstractVector{<:AbstractFloat}})

For an observer whose geocentric position and velocity are known,
prepare star-independent astrometry parameters for transformations
between ICRS and GCRS.  The Earth ephemeris is from ERFA models.

The parameters produced by this function are required in the space
motion, parallax, light deflection and aberration parts of the
astrometric transformation chain.

# Input

 - `date1`  -- TDB as a 2-part...
 - `date2`  -- ...Julian Date (Note 1)
 - `pv`     -- observer's geocentric pos/vel (Note 3)

# Output

 - 'astrom` -- star-independent astrometry parameters:

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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

2) All the vectors are with respect to BCRS axes.

3) The observer's position and velocity pv are geocentric but with
   respect to BCRS axes, and in units of m and m/s.  No assumptions
   are made about proximity to the Earth, and the function can be used
   for deep space applications as well as Earth orbit and terrestrial.

4) In cases where the caller wishes to supply his own Earth ephemeris,
   the function eraApcs can be used instead of the present function.

5) This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

		functions         observer        transformation

	 eraApcg eraApcg13    geocentric      ICRS <-> GCRS
	 eraApci eraApci13    terrestrial     ICRS <-> CIRS
	 eraApco eraApco13    terrestrial     ICRS <-> observed
	 eraApcs eraApcs13    space           ICRS <-> GCRS
	 eraAper eraAper13    terrestrial     update Earth rotation
	 eraApio eraApio13    terrestrial     CIRS <-> observed

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion, parallax,
   light deflection, and aberration.  From GCRS to CIRS comprises
   frame bias and precession-nutation.  From CIRS to observed takes
   account of Earth rotation, polar motion, diurnal aberration and
   parallax (unless subsumed into the ICRS <-> GCRS transformation),
   and atmospheric refraction.

6) The context structure astrom produced by this function is used by
   eraAtciq* and eraAticq*.
"""
function apcs13(day1::AbstractFloat, day2::AbstractFloat,
	pv::AbstractVector{<:AbstractVector{<:AbstractFloat}})
	#  Earth barycentric and heliocentric position & velocity (AU, AU/day).
	ehpv, ebpv = epv00(day1, day2)
	#  Compute the star-independent astrometry parameters.
	apcs(day1, day2, pv, ebpv, ehpv[1])
end

"""
	aper(θ::AbstractFloat, a::Astrom)

In the star-independent astrometry parameters, update only the Earth
rotation angle, supplied by the caller explicitly.

# Input

 - `theta`  -- Earth rotation angle (radians, Note 2)
 - `astrom` -- star-independent astrometry parameters:

# Output

 - `astrom` -- star-independent astrometry parameters:

# Note

1) This function exists to enable sidereal-tracking applications to
   avoid wasteful recomputation of the bulk of the astrometry
   parameters: only the Earth rotation is updated.

2) For targets expressed as equinox based positions, such as classical
   geocentric apparent (RA,Dec), the supplied theta can be Greenwich
   apparent sidereal time rather than Earth rotation angle.

3) The function eraAper13 can be used instead of the present function,
   and starts from UT1 rather than ERA itself.

4) This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

		functions         observer        transformation

	 eraApcg eraApcg13    geocentric      ICRS <-> GCRS
	 eraApci eraApci13    terrestrial     ICRS <-> CIRS
	 eraApco eraApco13    terrestrial     ICRS <-> observed
	 eraApcs eraApcs13    space           ICRS <-> GCRS
	 eraAper eraAper13    terrestrial     update Earth rotation
	 eraApio eraApio13    terrestrial     CIRS <-> observed

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion, parallax,
   light deflection, and aberration.  From GCRS to CIRS comprises
   frame bias and precession-nutation.  From CIRS to observed takes
   account of Earth rotation, polar motion, diurnal aberration and
   parallax (unless subsumed into the ICRS <-> GCRS transformation),
   and atmospheric refraction.

"""
function aper(θ::AbstractFloat, a::Astrom)
	Astrom(a.pmt, a.eb, a.eh, a.em, a.v, a.bm1, a.bpn, a.along, a.phi,
		a.xpl, a.ypl, a.sphi, a.cphi, a.diurab, θ+a.along, a.refa, a.refb)
end

"""
	aper13(day1::AbstractFloat, day2::AbstractFloat, a::Astrom)

In the star-independent astrometry parameters, update only the
Earth rotation angle.  The caller provides UT1, (n.b. not UTC).

# Input

 - `ut11`   -- UT1 as a 2-part...
 - `ut12`   -- ...Julian Date (Note 1)
 - `astrom` -- star-independent astrometry parameters:

# Output

 - `astrom` -- star-independent astrometry parameters:

# Note

1) The UT1 date (n.b. not UTC) ut11+ut12 is a Julian Date, apportioned
   in any convenient way between the arguments ut11 and ut12.  For
   example, JD(UT1)=2450123.7 could be expressed in any of these ways,
   among others:

		  ut11           ut12

	   2450123.7           0.0       (JD method)
	   2451545.0       -1421.3       (J2000 method)
	   2400000.5       50123.2       (MJD method)
	   2450123.5           0.2       (date & time method)

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 and MJD methods are good compromises between
   resolution and convenience.  The date & time method is best matched
   to the algorithm used: maximum precision is delivered when the ut11
   argument is for 0hrs UT1 on the day in question and the ut12
   argument lies in the range 0 to 1, or vice versa.

2) If the caller wishes to provide the Earth rotation angle itself,
   the function eraAper can be used instead.  One use of this
   technique is to substitute Greenwich apparent sidereal time and
   thereby to support equinox based transformations directly.

3) This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

		functions         observer        transformation

	 eraApcg eraApcg13    geocentric      ICRS <-> GCRS
	 eraApci eraApci13    terrestrial     ICRS <-> CIRS
	 eraApco eraApco13    terrestrial     ICRS <-> observed
	 eraApcs eraApcs13    space           ICRS <-> GCRS
	 eraAper eraAper13    terrestrial     update Earth rotation
	 eraApio eraApio13    terrestrial     CIRS <-> observed

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion, parallax,
   light deflection, and aberration.  From GCRS to CIRS comprises
   frame bias and precession-nutation.  From CIRS to observed takes
   account of Earth rotation, polar motion, diurnal aberration and
   parallax (unless subsumed into the ICRS <-> GCRS transformation),
   and atmospheric refraction.

"""
function aper13(day1::AbstractFloat, day2::AbstractFloat, a::Astrom)
	aper(era00(day1, day2), a)
end

"""
	apio(sp::AbstractFloat, θ::AbstractFloat, elong::AbstractFloat, ϕ::AbstractFloat, hm::AbstractFloat,
		 xp::AbstractFloat, yp::AbstractFloat, refa::AbstractFloat, refb::AbstractFloat, a::Astrom)

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between CIRS and observed coordinates.
The caller supplies the Earth orientation information and the
refraction constants as well as the site coordinates.

# Input

 - `sp`     -- the TIO locator s' (radians, Note 1)
 - `theta`  -- Earth rotation angle (radians)
 - `elong`  -- longitude (radians, east +ve, Note 2)
 - `phi`    -- geodetic latitude (radians, Note 2)
 - `hm`     -- height above ellipsoid (m, geodetic Note 2)
 - `xp,yp`  -- polar motion coordinates (radians, Note 3)
 - `refa`   -- refraction constant A (radians, Note 4)
 - `refb`   -- refraction constant B (radians, Note 4)

# Output

 - `astrom` -- star-independent astrometry parameters:

# Note

1) sp, the TIO locator s', is a tiny quantity needed only by the most
   precise applications.  It can either be set to zero or predicted
   using the ERFA function eraSp00.

2) The geographical coordinates are with respect to the ERFA_WGS84
   reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN: the
   longitude required by the present function is east-positive
   (i.e. right-handed), in accordance with geographical convention.

3) The polar motion xp,yp can be obtained from IERS bulletins.  The
   values are the coordinates (in radians) of the Celestial
   Intermediate Pole with respect to the International Terrestrial
   Reference System (see IERS Conventions 2003), measured along the
   meridians 0 and 90 deg west respectively.  For many applications,
   xp and yp can be set to zero.

   Internally, the polar motion is stored in a form rotated onto the
   local meridian.

4) The refraction constants refa and refb are for use in a dZ =
   A*tan(Z)+B*tan^3(Z) model, where Z is the observed (i.e. refracted)
   zenith distance and dZ is the amount of refraction.

5) It is advisable to take great care with units, as even unlikely
   values of the input parameters are accepted and processed in
   accordance with the models used.

6) In cases where the caller does not wish to provide the Earth
   rotation information and refraction constants, the function
   eraApio13 can be used instead of the present function.  This starts
   from UTC and weather readings etc. and computes suitable values
   using other ERFA functions.

7) This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

		functions         observer        transformation

	 eraApcg eraApcg13    geocentric      ICRS <-> GCRS
	 eraApci eraApci13    terrestrial     ICRS <-> CIRS
	 eraApco eraApco13    terrestrial     ICRS <-> observed
	 eraApcs eraApcs13    space           ICRS <-> GCRS
	 eraAper eraAper13    terrestrial     update Earth rotation
	 eraApio eraApio13    terrestrial     CIRS <-> observed

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion, parallax,
   light deflection, and aberration.  From GCRS to CIRS comprises
   frame bias and precession-nutation.  From CIRS to observed takes
   account of Earth rotation, polar motion, diurnal aberration and
   parallax (unless subsumed into the ICRS <-> GCRS transformation),
   and atmospheric refraction.

8) The context structure astrom produced by this function is used by
   eraAtioq and eraAtoiq.
"""
function apio(sp::AbstractFloat, θ::AbstractFloat, elong::AbstractFloat, ϕ::AbstractFloat, hm::AbstractFloat,
	xp::AbstractFloat, yp::AbstractFloat, refa::AbstractFloat, refb::AbstractFloat, a::Astrom)
	#  Form the rotation matrix, CIRS to apparent (HA, Dec).
	r = Rz(elong) * Rx(-yp) * Ry(-xp) * Rz(θ+sp)
	#  Solve for local Earth rotation angle.
	eral = r[1, 1] != 0.0 || r[1, 2] != 0.0 ? atan(r[1, 2], r[1, 1]) : 0.0
	#  Solve for polar motion (x, y) with respect to local meridian.
	xpl = atan(r[1, 3], norm(r[1, 1:2]))
	ypl = r[2, 3] != 0.0 || r[3, 3] != 0.0 ? -atan(r[2, 3], r[3, 3]) : 0.0
	#  Adjust longitude.
	along = anpm(eral - θ)
	#  Functions of latitude
	sphi, cphi = sincos(ϕ)
	#  Observer's geocentric position and velocity (m, m/s, CIRS).
	pv = pvtob(elong, ϕ, hm, xp, yp, sp, θ)
	#  Magnitude of diurnal aberration vector.
	diurab = norm(pv[2][1:2])/LIGHTSPEED
	Astrom(a.pmt, a.eb, a.eh, a.em, a.v, a.bm1, a.bpn, along, a.phi,
		xpl, ypl, sphi, cphi, diurab, eral, refa, refb)
end

"""
	 apio13(day1::AbstractFloat, day2::AbstractFloat, dut1::AbstractFloat, elong::AbstractFloat, ϕ::AbstractFloat,
			hm::AbstractFloat, xp::AbstractFloat, yp::AbstractFloat, phpa::AbstractFloat, tc::AbstractFloat,
			rh::AbstractFloat, wl::AbstractFloat, a::Astrom)

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between CIRS and observed coordinates.
The caller supplies UTC, site coordinates, ambient air conditions and
observing wavelength.

# Input

 - `utc1`   -- UTC as a 2-part...
 - `utc2`   -- ...quasi Julian Date (Notes 1,2)
 - `dut1`   -- UT1-UTC (seconds)
 - `elong`  -- longitude (radians, east +ve, Note 3)
 - `phi`    -- geodetic latitude (radians, Note 3)
 - `hm`     -- height above ellipsoid (m, geodetic Notes 4,6)
 - `xp,yp`  -- polar motion coordinates (radians, Note 5)
 - `phpa`   -- pressure at the observer (hPa = mB, Note 6)
 - `tc`     -- ambient temperature at the observer (deg C)
 - `rh`     -- relative humidity at the observer (range 0-1)
 - `wl`     -- wavelength (micrometers, Note 7)

# Output

 - `astrom` -- star-independent astrometry parameters:

# Note

1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
   convenient way between the two arguments, for example where utc1 is
   the Julian Day Number and utc2 is the fraction of a day.

   However, JD cannot unambiguously represent UTC during a leap second
   unless special measures are taken.  The convention in the present
   function is that the JD day represents UTC days whether the length
   is 86399, 86400 or 86401 SI seconds.

   Applications should use the function eraDtf2d to convert from
   calendar date and time of day into 2-part quasi Julian Date, as it
   implements the leap-second-ambiguity convention just described.

2) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale or that are too far in the future to
   be trusted.  See eraDat for further details.

3) UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
   one second at the end of each positive UTC leap second, introduced
   in order to keep UT1-UTC within +/- 0.9s.  n.b. This practice is
   under review, and in the future UT1-UTC may grow essentially
   without limit.

4) The geographical coordinates are with respect to the ERFA_WGS84
   reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN: the
   longitude required by the present function is east-positive
   (i.e. right-handed), in accordance with geographical convention.

5) The polar motion xp,yp can be obtained from IERS bulletins.  The
   values are the coordinates (in radians) of the Celestial
   Intermediate Pole with respect to the International Terrestrial
   Reference System (see IERS Conventions 2003), measured along the
   meridians 0 and 90 deg west respectively.  For many applications,
   xp and yp can be set to zero.

   Internally, the polar motion is stored in a form rotated onto the
   local meridian.

6) If hm, the height above the ellipsoid of the observing station in
   meters, is not known but phpa, the pressure in hPa (=mB), is
   available, an adequate estimate of hm can be obtained from the
   expression

		  hm = -29.3 * tsl * log ( phpa / 1013.25 );

   where tsl is the approximate sea-level air temperature in K (See
   Astrophysical Quantities, C.W.Allen, 3rd edition, section 52).
   Similarly, if the pressure phpa is not known, it can be estimated
   from the height of the observing station, hm, as follows:

		  phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

   Note, however, that the refraction is nearly proportional to the
   pressure and that an accurate phpa value is important for precise
   work.

7) The argument wl specifies the observing wavelength in micrometers.
   The transition from optical to radio is assumed to occur at 100
   micrometers (about 3000 GHz).

8) It is advisable to take great care with units, as even unlikely
   values of the input parameters are accepted and processed in
   accordance with the models used.

9) In cases where the caller wishes to supply his own Earth rotation
   information and refraction constants, the function eraApc can be
   used instead of the present function.

10) This is one of several functions that inserts into the astrom
	structure star-independent parameters needed for the chain of
	astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

	The various functions support different classes of observer and
	portions of the transformation chain:

		functions         observer        transformation

	 eraApcg eraApcg13    geocentric      ICRS <-> GCRS
	 eraApci eraApci13    terrestrial     ICRS <-> CIRS
	 eraApco eraApco13    terrestrial     ICRS <-> observed
	 eraApcs eraApcs13    space           ICRS <-> GCRS
	 eraAper eraAper13    terrestrial     update Earth rotation
	 eraApio eraApio13    terrestrial     CIRS <-> observed

	Those with names ending in "13" use contemporary ERFA models to
	compute the various ephemerides.  The others accept ephemerides
	supplied by the caller.

	The transformation from ICRS to GCRS covers space motion,
	parallax, light deflection, and aberration.  From GCRS to CIRS
	comprises frame bias and precession-nutation.  From CIRS to
	observed takes account of Earth rotation, polar motion, diurnal
	aberration and parallax (unless subsumed into the ICRS <-> GCRS
	transformation), and atmospheric refraction.

11) The context structure astrom produced by this function is used
	by eraAtioq and eraAtoiq.

"""
function apio13(day1::AbstractFloat, day2::AbstractFloat, dut1::AbstractFloat, elong::AbstractFloat,
	ϕ::AbstractFloat, hm::AbstractFloat, xp::AbstractFloat, yp::AbstractFloat, phpa::AbstractFloat,
	tc::AbstractFloat, rh::AbstractFloat, wl::AbstractFloat, a::Astrom)
	#  TIO locator s'.
	sp = sp00(taitt(utctai(day1, day2)...)...)
	#  Earth rotation angle.
	θ = era00(utcut1(day1, day2, dut1)...)
	#  Refraction constants A and B.
	refa, refb = refco(phpa, tc, rh, wl)
	#  CIRS <-> observed astrometry parameters.
	apio(sp, θ, elong, ϕ, hm, xp, yp, refa, refb, a)
end

"""
	atcc13(rc::AbstractFloat, dc::AbstractFloat, pr::AbstractFloat, pd::AbstractFloat, px::AbstractFloat,
		   rv::AbstractFloat, day1::AbstractFloat, day2::AbstractFloat)

Transform a star's ICRS catalog entry (epoch J2000.0) into ICRS
astrometric place.

# Input

 - `rc`    -- ICRS right ascension at J2000.0 (radians, Note 1)
 - `dc`    -- ICRS declination at J2000.0 (radians, Note 1)
 - `pr`    -- RA proper motion (radians/year, Note 2)
 - `pd`    -- Dec proper motion (radians/year)
 - `px`    -- parallax (arcsec)
 - `rv`    -- radial velocity (km/s, +ve if receding)
 - `date1` -- TDB as a 2-part...
 - `date2` -- ...Julian Date (Note 3)

# Output

 - ra, da` -- ICRS astrometric RA,Dec (radians)

# Note

1) Star data for an epoch other than J2000.0 (for example from the
   Hipparcos catalog, which has an epoch of J1991.25) will require a
   preliminary call to eraPmsafe before use.

2) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

3) The TDB date date1+date2 is a Julian Date, apportioned in any
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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

"""
function atcc13(rc::AbstractFloat, dc::AbstractFloat, pr::AbstractFloat, pd::AbstractFloat,
	px::AbstractFloat, rv::AbstractFloat, day1::AbstractFloat, day2::AbstractFloat)
	#  The transformation parameters and catalog ICRS (epoch j2000.0)
	#  to astrometric.
	atccq(rc, dc, pr, pd, px, rv, apci13(day1, day2)[1])
end

"""
	atccq(rc::AbstractFloat, dc::AbstractFloat, pr::AbstractFloat, pd::AbstractFloat, px::AbstractFloat,
		  rv::AbstractFloat, a::Astrom)

Quick transformation of a star's ICRS catalog entry (epoch J2000.0)
into ICRS astrometric place, given precomputed star-independent
astrometry parameters.

Use of this function is appropriate when efficiency is important and
where many star positions are to be transformed for one date.  The
star-independent parameters can be obtained by calling one of the
functions eraApci[13], eraApcg[13], eraApco[13] or eraApcs[13].

If the parallax and proper motions are zero the transformation has no
effect.

# Input

 - `rc,dc   -- ICRS RA,Dec at J2000.0 (radians)
 - `pr      -- RA proper motion (radians/year, Note 3)
 - `pd      -- Dec proper motion (radians/year)
 - `px      -- parallax (arcsec)
 - `rv      -- radial velocity (km/s, +ve if receding)
 - `astrom` -- star-independent astrometry parameters:

# Output

 - `ra, dec` -- RA, Dec of computed place

# Note

1) All the vectors are with respect to BCRS axes.

2) Star data for an epoch other than J2000.0 (for example from the
   Hipparcos catalog, which has an epoch of J1991.25) will require a
   preliminary call to eraPmsafe before use.

3) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

"""
function atccq(rc::AbstractFloat, dc::AbstractFloat, pr::AbstractFloat, pd::AbstractFloat, px::AbstractFloat,
	rv::AbstractFloat, a::Astrom)
	#  Proper motion and parallax, giving BCRS coordinate direction, and IRCS
	#  astrometric RA, Dec.
	ra, dec = c2s(pmpx(rc, dc, pr, pd, px, rv, a.pmt, a.eb))
	(anp(ra), dec)
end

"""
	atci13(rc::AbstractFloat, dc::AbstractFloat, pr::AbstractFloat, pd::AbstractFloat, px::AbstractFloat,
		   rv::AbstractFloat, day1::AbstractFloat, day2::AbstractFloat)

Transform ICRS star data, epoch J2000.0, to CIRS.

# Input

 - `rc`    -- ICRS right ascension at J2000.0 (radians, Note 1)
 - `dc`    -- ICRS declination at J2000.0 (radians, Note 1)
 - `pr`    -- RA proper motion (radians/year, Note 2)
 - `pd`    -- Dec proper motion (radians/year)
 - `px`    -- parallax (arcsec)
 - `rv`    -- radial velocity (km/s, +ve if receding)
 - `date1` -- TDB as a 2-part...
 - `date2` -- ...Julian Date (Note 3)

# Output

  - `ri, di` -- CIRS geocentric RA,Dec (radians)
  - `eo`     -- equation of the origins (ERA-GST, Note 5)

# Note

1) Star data for an epoch other than J2000.0 (for example from the
   Hipparcos catalog, which has an epoch of J1991.25) will require a
   preliminary call to eraPmsafe before use.

2) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

3) The TDB date date1+date2 is a Julian Date, apportioned in any
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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

4) The available accuracy is better than 1 milliarcsecond, limited
   mainly by the precession-nutation model that is used, namely IAU
   2000A/2006.  Very close to solar system bodies, additional errors
   of up to several milliarcseconds can occur because of unmodeled
   light deflection; however, the Sun's contribution is taken into
   account, to first order.  The accuracy limitations of the ERFA
   function eraEpv00 (used to compute Earth position and velocity) can
   contribute aberration errors of up to 5 microarcseconds.  Light
   deflection at the Sun's limb is uncertain at the 0.4 mas level.

5) Should the transformation to (equinox based) apparent place be
   required rather than (CIO based) intermediate place, subtract the
   equation of the origins from the returned right ascension: RA = RI
   - EO. (The eraAnp function can then be applied, as required, to
   keep the result in the conventional 0-2pi range.)
"""
function atci13(rc::AbstractFloat, dc::AbstractFloat, pr::AbstractFloat, pd::AbstractFloat,
	px::AbstractFloat, rv::AbstractFloat, day1::AbstractFloat, day2::AbstractFloat)
	#  The transformation parameters.
	a, eo = apci13(day1, day2)
	#  ICRS (epoch J2000.0) to CIRS.
	ra, dec = values(atciq(rc, dc, pr, pd, px, rv, a))
	(ra = ra, dec = dec, eo = eo)
end

"""
	atciq(rc::AbstractFloat, dc::AbstractFloat, pr::AbstractFloat, pd::AbstractFloat, px::AbstractFloat,
		  rv::AbstractFloat, a::Astrom)

Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
star-independent astrometry parameters.

Use of this function is appropriate when efficiency is important and
where many star positions are to be transformed for one date.  The
star-independent parameters can be obtained by calling one of the
functions eraApci[13], eraApcg[13], eraApco[13] or eraApcs[13].

If the parallax and proper motions are zero the eraAtciqz function can
be used instead.

# Input

 - `rc, dc` -- ICRS RA,Dec at J2000.0 (radians, Note 1)
 - `pr`     -- RA proper motion (radians/year, Note 2)
 - `pd`     -- Dec proper motion (radians/year)
 - `px`     -- parallax (arcsec)
 - `rv`     -- radial velocity (km/s, +ve if receding)
 - `astrom` -- star-independent astrometry parameters:

# Output

 - `ra, dec` -- RA, Dec of transformed coordinates.

# Note

1) Star data for an epoch other than J2000.0 (for example from the
   Hipparcos catalog, which has an epoch of J1991.25) will require a
   preliminary call to eraPmsafe before use.

2) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

"""
function atciq(rc::AbstractFloat, dc::AbstractFloat, pr::AbstractFloat, pd::AbstractFloat,
	px::AbstractFloat, rv::AbstractFloat, a::Astrom)
	#  Proper motion and parallax, giving BCRS coordinate direction.
	#  Light deflection by the Sun, giving BCRS natural direction.
	#  Aberration, giving GCRS proper direction.
	#  Bias-precession-Nutation, giving CIRS proper direction.
	ri, di = c2s(a.bpn*ab(ldsun(pmpx(rc, dc, pr, pd, px, rv, a.pmt, a.eb),
			a.eh, a.em), a.v, a.em, a.bm1))
	(ra = mod2pi(ri), dec = di)
end

"""
	atciqn(rc::AbstractFloat, dc::AbstractFloat, pr::AbstractFloat, pd::AbstractFloat, px::AbstractFloat,
		   rv::AbstractFloat, a::Astrom, n::Int, b::AbstractVector{Ldbody})

Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
star-independent astrometry parameters plus a list of light-
deflecting bodies.

Use of this function is appropriate when efficiency is important and
where many star positions are to be transformed for one date.  The
star-independent parameters can be obtained by calling one of the
functions eraApci[13], eraApcg[13], eraApco[13] or eraApcs[13].

If the only light-deflecting body to be taken into account is the Sun,
the eraAtciq function can be used instead.  If in addition the
parallax and proper motions are zero, the eraAtciqz function can be
used.

# Input

 - `rc, dc` -- ICRS RA,Dec at J2000.0 (radians)
 - `pr`     -- RA proper motion (radians/year, Note 3)
 - `pd`     -- Dec proper motion (radians/year)
 - `px`     -- parallax (arcsec)
 - `rv`     -- radial velocity (km/s, +ve if receding)
 - `astrom` -- star-independent astrometry parameters:

# Output

 - `ri, di` -- CIRS RA,Dec (radians)

# Note

1) Star data for an epoch other than J2000.0 (for example from the
   Hipparcos catalog, which has an epoch of J1991.25) will require a
   preliminary call to eraPmsafe before use.

2) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

3) The struct b contains n entries, one for each body to be
   considered.  If n = 0, no gravitational light deflection will be
   applied, not even for the Sun.

4) The struct b should include an entry for the Sun as well as for any
   planet or other body to be taken into account.  The entries should
   be in the order in which the light passes the body.

5) In the entry in the b struct for body i, the mass parameter b[i].bm
   can, as required, be adjusted in order to allow for such effects as
   quadrupole field.

6) The deflection limiter parameter b[i].dl is phi^2/2, where phi is
   the angular separation (in radians) between star and body at which
   limiting is applied.  As phi shrinks below the chosen threshold,
   the deflection is artificially reduced, reaching zero for phi = 0.
   Example values suitable for a terrestrial observer, together with
   masses, are as follows:

	  body i     b[i].bm        b[i].dl

	  Sun        1.0            6e-6
	  Jupiter    0.00095435     3e-9
	  Saturn     0.00028574     3e-10

7) For efficiency, validation of the contents of the b array is
   omitted.  The supplied masses must be greater than zero, the
   position and velocity vectors must be right, and the deflection
   limiter greater than zero.
"""
function atciqn(rc::AbstractFloat, dc::AbstractFloat, pr::AbstractFloat, pd::AbstractFloat,
	px::AbstractFloat, rv::AbstractFloat, a::Astrom, n::Int,
	b::AbstractVector{Ldbody})
	#  Proper motion and parallax, giving BCRS coordinate direction.
	#  Light deflection, giving natural direction.
	#  Aberration, giving GCRS proper direction.
	#  Bias-precession-nutation, giving CIRS proper direction.
	#  CIRS RA, Dec
	ra, dec = c2s(a.bpn*ab(
		ldn(n, b, a.eb, pmpx(rc, dc, pr, pd, px, rv, a.pmt, a.eb)),
		a.v, a.em, a.bm1))
	(ra = mod2pi(ra), dec = dec)
end

"""
	atciqz(rc::AbstractFloat, dc::AbstractFloat, a::Astrom)

Quick ICRS to CIRS transformation, given precomputed star- independent
astrometry parameters, and assuming zero parallax and proper motion.

Use of this function is appropriate when efficiency is important and
where many star positions are to be transformed for one date.  The
star-independent parameters can be obtained by calling one of the
functions eraApci[13], eraApcg[13], eraApco[13] or eraApcs[13].

The corresponding function for the case of non-zero parallax and
proper motion is eraAtciq.

# Input

 - `rc`    -- ICRS astrometric RA (radians)
 - `dc`    -- ICRS astrometric Dec (radians)
 - `astrom` -- star-independent astrometry parameters:

# Output

 - `ri, di` -- CIRS RA,Dec (radians)

# Note

   All the vectors are with respect to BCRS axes.

# References

Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to the
Astronomical Almanac, 3rd ed., University Science Books (2013).

Klioner, Sergei A., "A practical relativistic model for micro-
arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).

"""
function atciqz(rc::AbstractFloat, dc::AbstractFloat, a::Astrom)
	#  BCRS coordinate direction (unit vector).
	#  Light deflection by the Sun, giving BCRS natural direction.
	#  Aberration, giving GCRS proper direction.
	#  Bias-precession-nutation, giving CIRS proper direction.
	#  CIRS to RA, Dec.
	ra, dec = c2s(a.bpn*ab(ldsun(s2c(rc, dc), a.eh, a.em), a.v, a.em, a.bm1))
	(ra = mod2pi(ra), dec = dec)
end

"""
	atco13(rc::AbstractFloat, dc::AbstractFloat, pr::AbstractFloat, pd::AbstractFloat, px::AbstractFloat,
		   rv::AbstractFloat, utc1::AbstractFloat, utc2::AbstractFloat, dut1::AbstractFloat,
		   elong::AbstractFloat, ϕ::AbstractFloat, rm::AbstractFloat, xp::AbstractFloat, yp::AbstractFloat,
		   phpa::AbstractFloat, tc::AbstractFloat, rh::AbstractFloat, wl::AbstractFloat)

ICRS RA,Dec to observed place.  The caller supplies UTC, site
coordinates, ambient air conditions and observing wavelength.

SOFA models are used for the Earth ephemeris, bias-precession-
nutation, Earth orientation and refraction.

# Input

 - `rc, dc` -- ICRS right ascension at J2000.0 (radians, Note 1)
 - `pr`     -- RA proper motion (radians/year, Note 2)
 - `pd`     -- Dec proper motion (radians/year)
 - `px`     -- parallax (arcsec)
 - `rv`     -- radial velocity (km/s, +ve if receding)
 - `utc1`   -- UTC as a 2-part...
 - `utc2`   -- ...quasi Julian Date (Notes 3-4)
 - `dut1`   -- UT1-UTC (seconds, Note 5)
 - `elong`  -- longitude (radians, east +ve, Note 6)
 - `phi`    -- latitude (geodetic, radians, Note 6)
 - `hm`     -- height above ellipsoid (m, geodetic, Notes 6,8)
 - `xp, yp` -- polar motion coordinates (radians, Note 7)
 - `phpa`   -- pressure at the observer (hPa = mB, Note 8)
 - `tc`     -- ambient temperature at the observer (deg C)
 - `rh`     -- relative humidity at the observer (range 0-1)
 - `wl`     -- wavelength (micrometers, Note 9)

# Output

 - `aob`    -- observed azimuth (radians: N=0,E=90)
 - `zob`    -- observed zenith distance (radians)
 - `hob`    -- observed hour angle (radians)
 - `dob`    -- observed declination (radians)
 - `rob`    -- observed right ascension (CIO-based, radians)
 - `eo`     -- equation of the origins (ERA-GST)

# Note

1) Star data for an epoch other than J2000.0 (for example from the
   Hipparcos catalog, which has an epoch of J1991.25) will require a
   preliminary call to eraPmsafe before use.

2) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

3) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
   convenient way between the two arguments, for example where utc1 is
   the Julian Day Number and utc2 is the fraction of a day.

   However, JD cannot unambiguously represent UTC during a leap second
   unless special measures are taken.  The convention in the present
   function is that the JD day represents UTC days whether the length
   is 86399, 86400 or 86401 SI seconds.

   Applications should use the function eraDtf2d to convert from
   calendar date and time of day into 2-part quasi Julian Date, as it
   implements the leap-second-ambiguity convention just described.

4) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale or that are too far in the future to
   be trusted.  See eraDat for further details.

5) UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
   one second at the end of each positive UTC leap second, introduced
   in order to keep UT1-UTC within +/- 0.9s.  n.b. This practice is
   under review, and in the future UT1-UTC may grow essentially
   without limit.

6) The geographical coordinates are with respect to the ERFA_WGS84
   reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN: the
   longitude required by the present function is east-positive
   (i.e. right-handed), in accordance with geographical convention.
   7) The polar motion xp,yp can be obtained from IERS bulletins.  The
   values are the coordinates (in radians) of the Celestial
   Intermediate Pole with respect to the International Terrestrial
   Reference System (see IERS Conventions 2003), measured along the
   meridians 0 and 90 deg west respectively.  For many applications,
   xp and yp can be set to zero.

8) If hm, the height above the ellipsoid of the observing station in
   meters, is not known but phpa, the pressure in hPa (=mB), is
   available, an adequate estimate of hm can be obtained from the
   expression

		  hm = -29.3 * tsl * log ( phpa / 1013.25 );

   where tsl is the approximate sea-level air temperature in K (See
   Astrophysical Quantities, C.W.Allen, 3rd edition, section 52).
   Similarly, if the pressure phpa is not known, it can be estimated
   from the height of the observing station, hm, as follows:

		  phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

   Note, however, that the refraction is nearly proportional to the
   pressure and that an accurate phpa value is important for precise
   work.

9) The argument wl specifies the observing wavelength in micrometers.
   The transition from optical to radio is assumed to occur at 100
   micrometers (about 3000 GHz).

10) The accuracy of the result is limited by the corrections for
	refraction, which use a simple A*tan(z) + B*tan^3(z) model.
	Providing the meteorological parameters are known accurately and
	there are no gross local effects, the predicted observed
	coordinates should be within 0.05 arcsec (optical) or 1 arcsec
	(radio) for a zenith distance of less than 70 degrees, better than
	30 arcsec (optical or radio) at 85 degrees and better than 20
	arcmin (optical) or 30 arcmin (radio) at the horizon.

	Without refraction, the complementary functions eraAtco13 and
	eraAtoc13 are self-consistent to better than 1 microarcsecond all
	over the celestial sphere.  With refraction included, consistency
	falls off at high zenith distances, but is still better than 0.05
	arcsec at 85 degrees.

11) "Observed" Az,ZD means the position that would be seen by a
	perfect geodetically aligned theodolite.  (Zenith distance is used
	rather than altitude in order to reflect the fact that no
	allowance is made for depression of the horizon.)  This is related
	to the observed HA,Dec via the standard rotation, using the
	geodetic latitude (corrected for polar motion), while the observed
	HA and RA are related simply through the Earth rotation angle and
	the site longitude.  "Observed" RA,Dec or HA,Dec thus means the
	position that would be seen by a perfect equatorial with its polar
	axis aligned to the Earth's axis of rotation.

12) It is advisable to take great care with units, as even unlikely
	values of the input parameters are accepted and processed in
	accordance with the models used.

"""
function atco13(rc::AbstractFloat, dc::AbstractFloat, pr::AbstractFloat, pd::AbstractFloat,
	px::AbstractFloat, rv::AbstractFloat, utc1::AbstractFloat, utc2::AbstractFloat,
	dut1::AbstractFloat, elong::AbstractFloat, ϕ::AbstractFloat, rm::AbstractFloat,
	xp::AbstractFloat, yp::AbstractFloat, phpa::AbstractFloat, tc::AbstractFloat,
	rh::AbstractFloat, wl::AbstractFloat)
	#  Star-independent astrometry parameters.
	a, eo = apco13(utc1, utc2, dut1, elong, ϕ, rm, xp, yp, phpa, tc, rh, wl)
	#  Transform ICRS to CIRS and CIRS to observed.
	azi, zen, ha, dec, ra = atioq(atciq(rc, dc, pr, pd, px, rv, a)..., a)
	(azi = azi, zen = zen, ha = ha, dec = dec, ra = ra, eo = eo)
end

"""
	atic13(ri::AbstractFloat, di::AbstractFloat, day1::AbstractFloat, day2::AbstractFloat)

Transform star RA,Dec from geocentric CIRS to ICRS astrometric.

# Input

 - `ri, di` -- CIRS geocentric RA,Dec (radians)
 - `date1`  -- TDB as a 2-part...
 - `date2`  -- ...Julian Date (Note 1)

# Output

 - `rc, dc` -- ICRS astrometric RA,Dec (radians)
 - `eo`     -- equation of the origins (ERA-GST, Note 4)

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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

2) Iterative techniques are used for the aberration and light
   deflection corrections so that the functions eraAtic13 (or
   eraAticq) and eraAtci13 (or eraAtciq) are accurate inverses; even
   at the edge of the Sun's disk the discrepancy is only about 1
   nanoarcsecond.

3) The available accuracy is better than 1 milliarcsecond, limited
   mainly by the precession-nutation model that is used, namely IAU
   2000A/2006.  Very close to solar system bodies, additional errors
   of up to several milliarcseconds can occur because of unmodeled
   light deflection; however, the Sun's contribution is taken into
   account, to first order.  The accuracy limitations of the ERFA
   function eraEpv00 (used to compute Earth position and velocity) can
   contribute aberration errors of up to 5 microarcseconds.  Light
   deflection at the Sun's limb is uncertain at the 0.4 mas level.

4) Should the transformation to (equinox based) J2000.0 mean place be
   required rather than (CIO based) ICRS coordinates, subtract the
   equation of the origins from the returned right ascension: RA = RI
   - EO.  (The eraAnp function can then be applied, as required, to
   keep the result in the conventional 0-2pi range.)

"""
function atic13(ri::AbstractFloat, di::AbstractFloat, day1::AbstractFloat, day2::AbstractFloat)
	a, eo = apci13(day1, day2)
	ra, dec = aticq(ri, di, a)
	(ra = ra, dec = dec, eo)
end

"""
	atic13(ri::AbstractFloat, di::AbstractFloat, day1::AbstractFloat, day2::AbstractFloat)

Quick CIRS RA,Dec to ICRS astrometric place, given the star-
independent astrometry parameters.

Use of this function is appropriate when efficiency is important and
where many star positions are all to be transformed for one date.  The
star-independent astrometry parameters can be obtained by calling one
of the functions eraApci[13], eraApcg[13], eraApco[13] or eraApcs[13].

# Input

 - `ri, di` -- CIRS RA,Dec (radians)
 - `astrom` -- star-independent astrometry parameters:

# Output

 - `rc, dc` -- ICRS astrometric RA,Dec (radians)

# Note

1) Only the Sun is taken into account in the light deflection
   correction.

2) Iterative techniques are used for the aberration and light
   deflection corrections so that the functions eraAtic13 (or
   eraAticq) and eraAtci13 (or eraAtciq) are accurate inverses; even
   at the edge of the Sun's disk the discrepancy is only about 1
   nanoarcsecond.
"""
function aticq(ri::AbstractFloat, di::AbstractFloat, a::Astrom)
	#  CIRS RA, Dec to cartesion.
	#  Bias-precession-nutation, giving GCRS proper direction.
	ppr = a.bpn'*s2c(ri, di)
	#  Aberration, giving GCRS natural direction.
	d, pnat, pco = zeros(Float64, 3), zeros(Float64, 3), zeros(Float64, 3)
	for j ∈ 1:3
		w = ppr .- d
		bf = copy(w) ./ norm(w)
		af = ab(bf, a.v, a.em, a.bm1)
		d = af .- bf
		w = ppr .- d
		pnat .= copy(w) ./ norm(w)
	end
	#  Light deflection by the Sun, giving BCRS coordinate direction.
	d = zeros(Float64, 3)
	for j ∈ 1:3
		w = pnat .- d
		bf = copy(w) ./ norm(w)
		af = ldsun(bf, a.eh, a.em)
		d = af .- bf
		w = pnat .- d
		pco .= copy(w) ./ norm(w)
	end
	#  ICRS astrometric RA, Dec.
	ra, dec = c2s(pco)
	(ra = mod2pi(ra), dec = dec)
end

"""
	aticqn(ri::AbstractFloat, di::AbstractFloat, a::Astrom, n::Int, b::AbstractVector{Ldbody})

Quick CIRS to ICRS astrometric place transformation, given the star-
independent astrometry parameters plus a list of light-deflecting
bodies.

Use of this function is appropriate when efficiency is important and
where many star positions are all to be transformed for one date.  The
star-independent astrometry parameters can be obtained by calling one
of the functions eraApci[13], eraApcg[13], eraApco[13] or eraApcs[13].

If the only light-deflecting body to be taken into account is the Sun,
the eraAticq function can be used instead.

# Input

 - `ri, di` -- CIRS RA,Dec (radians)
 - `astrom` -- star-independent astrometry parameters:
 - `n`      -- number of bodies (Note 3)
 - `b`      -- data for each of the n bodies (Notes 3,4):

# Output

 - `rc, dc` -- ICRS astrometric RA,Dec (radians)

# Note

1) Iterative techniques are used for the aberration and light
   deflection corrections so that the functions eraAticqn and
   eraAtciqn are accurate inverses; even at the edge of the Sun's disk
   the discrepancy is only about 1 nanoarcsecond.

2) If the only light-deflecting body to be taken into account is the
   Sun, the eraAticq function can be used instead.

3) The struct b contains n entries, one for each body to be
   considered.  If n = 0, no gravitational light deflection will be
   applied, not even for the Sun.

4) The struct b should include an entry for the Sun as well as for any
   planet or other body to be taken into account.  The entries should
   be in the order in which the light passes the body.

5) In the entry in the b struct for body i, the mass parameter b[i].bm
   can, as required, be adjusted in order to allow for such effects as
   quadrupole field.

6) The deflection limiter parameter b[i].dl is phi^2/2, where phi is
   the angular separation (in radians) between star and body at which
   limiting is applied.  As phi shrinks below the chosen threshold,
   the deflection is artificially reduced, reaching zero for phi = 0.
   Example values suitable for a terrestrial observer, together with
   masses, are as follows:

	  body i     b[i].bm        b[i].dl

	  Sun        1.0            6e-6
	  Jupiter    0.00095435     3e-9
	  Saturn     0.00028574     3e-10

7) For efficiency, validation of the contents of the b array is
   omitted.  The supplied masses must be greater than zero, the
   position and velocity vectors must be right, and the deflection
   limiter greater than zero.
"""
function aticqn(ri::AbstractFloat, di::AbstractFloat, a::Astrom, n::Int,
      b::AbstractVector{Ldbody})
	#  CIRS RA, Dec to cartesian.
	#  Bias-precession-nutation, giving GCRS proper direction.
	ppr = a.bpn'*s2c(ri, di)
	#  Aberration, giving GCRS natural direction.
	d, pnat, pco = zeros(Float64, 3), zeros(Float64, 3), zeros(Float64, 3)
	for j ∈ 1:2
		bf = (ppr .- d) ./ norm(ppr .- d)
		af = ab(bf, a.v, a.em, a.bm1)
		d .= af .- bf
		pnat .= (ppr .- d) ./ norm(ppr .- d)
	end
	#  Light deflection, giving BCRS coordinate direction.
	d = zeros(Float64, 3)
	for j ∈ 1:5
		bf = (pnat .- d) ./ norm(pnat .- d)
		af = ldn(n, b, a.eb, bf)
		d .= af .- bf
		pco .= (pnat .- d) ./ norm(pnat .- d)
	end
	rc, dc = c2s(pco)
	(ra = mod2pi(rc), dec = dc)
end

"""
	atio13(ri::AbstractFloat, di::AbstractFloat, utc1::AbstractFloat, utc2::AbstractFloat, dut1::AbstractFloat,
		   elong::AbstractFloat, ϕ::AbstractFloat, hm::AbstractFloat, xp::AbstractFloat, yp::AbstractFloat,
		   phpa::AbstractFloat, tc::AbstractFloat, rh::AbstractFloat, wl::AbstractFloat)

CIRS RA,Dec to observed place.  The caller supplies UTC, site
coordinates, ambient air conditions and observing wavelength.

# Input

 - `ri`     -- CIRS right ascension (CIO-based, radians)
 - `di`     -- CIRS declination (radians)
 - `utc1`   -- UTC as a 2-part...
 - `utc2`   -- ...quasi Julian Date (Notes 1,2)
 - `dut1`   -- UT1-UTC (seconds, Note 3)
 - `elong`  -- longitude (radians, east +ve, Note 4)
 - `phi`    -- geodetic latitude (radians, Note 4)
 - `hm`     -- height above ellipsoid (m, geodetic Notes 4,6)
 - `xp, yp` -- polar motion coordinates (radians, Note 5)
 - `phpa`   -- pressure at the observer (hPa = mB, Note 6)
 - `tc`     -- ambient temperature at the observer (deg C)
 - `rh`     -- relative humidity at the observer (range 0-1)
 - `wl`     -- wavelength (micrometers, Note 7)

# Output

 - `aob`    -- observed azimuth (radians: N=0,E=90)
 - `zob`    -- observed zenith distance (radians)
 - `hob`    -- observed hour angle (radians)
 - `dob`    -- observed declination (radians)
 - `rob`    -- observed right ascension (CIO-based, radians)

# Note

1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
   convenient way between the two arguments, for example where utc1 is
   the Julian Day Number and utc2 is the fraction of a day.

   However, JD cannot unambiguously represent UTC during a leap second
   unless special measures are taken.  The convention in the present
   function is that the JD day represents UTC days whether the length
   is 86399, 86400 or 86401 SI seconds.

   Applications should use the function eraDtf2d to convert from
   calendar date and time of day into 2-part quasi Julian Date, as it
   implements the leap-second-ambiguity convention just described.

2) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale or that are too far in the future to
   be trusted.  See eraDat for further details.

3) UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
   one second at the end of each positive UTC leap second, introduced
   in order to keep UT1-UTC within +/- 0.9s.  n.b. This practice is
   under review, and in the future UT1-UTC may grow essentially
   without limit.

4) The geographical coordinates are with respect to the ERFA_WGS84
   reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN: the
   longitude required by the present function is east-positive
   (i.e. right-handed), in accordance with geographical convention.

5) The polar motion xp,yp can be obtained from IERS bulletins.  The
   values are the coordinates (in radians) of the Celestial
   Intermediate Pole with respect to the International Terrestrial
   Reference System (see IERS Conventions 2003), measured along the
   meridians 0 and 90 deg west respectively.  For many applications,
   xp and yp can be set to zero.

6) If hm, the height above the ellipsoid of the observing station in
   meters, is not known but phpa, the pressure in hPa (=mB), is
   available, an adequate estimate of hm can be obtained from the
   expression

		  hm = -29.3 * tsl * log ( phpa / 1013.25 );

   where tsl is the approximate sea-level air temperature in K (See
   Astrophysical Quantities, C.W.Allen, 3rd edition, section 52).
   Similarly, if the pressure phpa is not known, it can be estimated
   from the height of the observing station, hm, as follows:

		  phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

   Note, however, that the refraction is nearly proportional to the
   pressure and that an accurate phpa value is important for precise
   work.

7) The argument wl specifies the observing wavelength in micrometers.
   The transition from optical to radio is assumed to occur at 100
   micrometers (about 3000 GHz).

8) "Observed" Az,ZD means the position that would be seen by a perfect
   geodetically aligned theodolite.  (Zenith distance is used rather
   than altitude in order to reflect the fact that no allowance is
   made for depression of the horizon.)  This is related to the
   observed HA,Dec via the standard rotation, using the geodetic
   latitude (corrected for polar motion), while the observed HA and RA
   are related simply through the Earth rotation angle and the site
   longitude.  "Observed" RA,Dec or HA,Dec thus means the position
   that would be seen by a perfect equatorial with its polar axis
   aligned to the Earth's axis of rotation.

9) The accuracy of the result is limited by the corrections for
   refraction, which use a simple A*tan(z) + B*tan^3(z) model.
   Providing the meteorological parameters are known accurately and
   there are no gross local effects, the predicted astrometric
   coordinates should be within 0.05 arcsec (optical) or 1 arcsec
   (radio) for a zenith distance of less than 70 degrees, better than
   30 arcsec (optical or radio) at 85 degrees and better than 20
   arcmin (optical) or 30 arcmin (radio) at the horizon.

10) The complementary functions eraAtio13 and eraAtoi13 are self-
	consistent to better than 1 microarcsecond all over the celestial
	sphere.

11) It is advisable to take great care with units, as even unlikely
	values of the input parameters are accepted and processed in
	accordance with the models used.
"""
function atio13(ri::AbstractFloat, di::AbstractFloat, utc1::AbstractFloat, utc2::AbstractFloat, dut1::AbstractFloat,
	elong::AbstractFloat, ϕ::AbstractFloat, hm::AbstractFloat, xp::AbstractFloat, yp::AbstractFloat,
	phpa::AbstractFloat, tc::AbstractFloat, rh::AbstractFloat, wl::AbstractFloat)
	#  Star-independent astrometry parameters for CIRS->observed.
	#  Transform CIRS to observed.
	aob, zob, hob, dob, rob = atioq(ri, di, apio13(
		utc1, utc2, dut1, elong, ϕ, hm, xp, yp, phpa, tc, rh, wl, Astrom()))
	(azi = aob, zen = zob, ha = hob, dec = dob, ra = rob)
end

"""
	atioq(ri::AbstractFloat, di::AbstractFloat, a::Astrom)

Use of this function is appropriate when efficiency is important and
where many star positions are all to be transformed for one date.  The
star-independent astrometry parameters can be obtained by calling
eraApio[13] or eraApco[13].

# Input

 - `ri`     -- CIRS right ascension
 - `di`     -- CIRS declination
 - `astrom` -- star-independent astrometry parameters:

# Output

 - `aob`    -- observed azimuth (radians: N=0,E=90)
 - `zob`    -- observed zenith distance (radians)
 - `hob`    -- observed hour angle (radians)
 - `dob`    -- observed declination (radians)
 - `rob`    -- observed right ascension (CIO-based, radians)

# Note

1) This function returns zenith distance rather than altitude in order
   to reflect the fact that no allowance is made for depression of the
   horizon.

2) The accuracy of the result is limited by the corrections for
   refraction, which use a simple A*tan(z) + B*tan^3(z) model.
   Providing the meteorological parameters are known accurately and
   there are no gross local effects, the predicted observed
   coordinates should be within 0.05 arcsec (optical) or 1 arcsec
   (radio) for a zenith distance of less than 70 degrees, better than
   30 arcsec (optical or radio) at 85 degrees and better than 20
   arcmin (optical) or 30 arcmin (radio) at the horizon.

   Without refraction, the complementary functions eraAtioq and
   eraAtoiq are self-consistent to better than 1 microarcsecond all
   over the celestial sphere.  With refraction included, consistency
   falls off at high zenith distances, but is still better than 0.05
   arcsec at 85 degrees.

3) It is advisable to take great care with units, as even unlikely
   values of the input parameters are accepted and processed in
   accordance with the models used.

4) The CIRS RA,Dec is obtained from a star catalog mean place by
   allowing for space motion, parallax, the Sun's gravitational lens
   effect, annual aberration and precession-nutation.  For star
   positions in the ICRS, these effects can be applied by means of the
   eraAtci13 (etc.) functions.  Starting from classical "mean place"
   systems, additional transformations will be needed first.

5) "Observed" Az,El means the position that would be seen by a perfect
   geodetically aligned theodolite.  This is obtained from the CIRS
   RA,Dec by allowing for Earth orientation and diurnal aberration,
   rotating from equator to horizon coordinates, and then adjusting
   for refraction.  The HA,Dec is obtained by rotating back into
   equatorial coordinates, and is the position that would be seen by a
   perfect equatorial with its polar axis aligned to the Earth's axis
   of rotation.  Finally, the RA is obtained by subtracting the HA
   from the local ERA.

6) The star-independent CIRS-to-observed-place parameters in ASTROM
   may be computed with eraApio[13] or eraApco[13].  If nothing has
   changed significantly except the time, eraAper[13] may be used to
   perform the requisite adjustment to the astrom structure.
"""
function atioq(ri::AbstractFloat, di::AbstractFloat, a::Astrom)
	#  Minimum cos(alt) and sin(alt) for refraction.
	CELMIN, SELMIN = 1e-6, 0.05
	#  CIRS Ra, Dec to cartesian and polar motion.
	hd = SMatrix{3, 3}([             cos(a.xpl) 0.0 sin(a.xpl);
		sin(a.xpl)*sin(a.ypl) cos(a.ypl) -cos(a.xpl)*sin(a.ypl);
		-sin(a.xpl)*cos(a.ypl) sin(a.ypl) cos(a.xpl)*cos(a.ypl)])*s2c(ri - a.eral, di)
	#  Diurnal aberration.
	hdt = (1.0 - a.diurab*hd[2]) .* (hd .+ SVector(0.0, a.diurab, 0.0))
	#  Cartesian -HA, Dec to cartesion Az, El (S=0, E=90).
	aet = SMatrix{3, 3}([a.sphi 0.0 -a.cphi; 0.0 1.0 0.0; a.cphi 0.0 a.sphi])*hdt
	#  Azimuth (N=0, E=90)
	azob = aet[1] != 0.0 || aet[2] != 0.0 ? atan(aet[2], -aet[1]) : 0.0
	####    Refraction    ####
	#  Cosine and sine of altitude, with precautions.
	r, z = maximum([norm(aet[1:2]) CELMIN; aet[3] SELMIN], dims = 2)
	#  A*tan(z) + B*tan^3(z) model, with Newton-Raphson correction.
	w = a.refb*(r/z)^2
	del = (a.refa + w)*(r/z)/(1.0 + (a.refa + 3*w)/z^2)
	#  Apply the change, giving observed vector.
	cosdel = 1.0 - del^2/2.0
	aeo = SVector((cosdel - del*z/r), cosdel - del*z/r, cosdel) .* aet .+ SVector(0.0, 0.0, del*r)
	#  Observed ZD.
	zdob = atan(norm(aeo[1:2]), aeo[3])
	#  Az/El vector to HA, Dec vector (both right-handed) and to spherical -HA, Dec.
	hmob, dcob = c2s(SMatrix{3, 3}([a.sphi 0.0 a.cphi; 0.0 1.0 0.0; -a.cphi 0.0 a.sphi])*aeo)
	#  Right ascension (with respect to CIO).
	raob = a.eral + hmob
	(azi = mod2pi(azob), zen = zdob, ha = -hmob, dec = dcob, ra = mod2pi(raob))
end

"""
   atoc13(tp::Char, ob1::AbstractFloat, ob2::AbstractFloat, utc1::AbstractFloat, utc2::AbstractFloat,
		  dut1::AbstractFloat, elong::AbstractFloat, ϕ::AbstractFloat, hm::AbstractFloat, xp::AbstractFloat,
		  yp::AbstractFloat, phpa::AbstractFloat, tc::AbstractFloat, rh::AbstractFloat, wl::AbstractFloat)

Observed place at a groundbased site to to ICRS astrometric RA,Dec.
The caller supplies UTC, site coordinates, ambient air conditions and
observing wavelength.

# Input

 - `type`   -- type of coordinates - "R", "H" or "A" (Notes 1,2)
 - `ob1`    -- observed Az, HA or RA (radians; Az is N=0,E=90)
 - `ob2`    -- observed ZD or Dec (radians)
 - `utc1`   -- UTC as a 2-part...
 - `utc2`   -- ...quasi Julian Date (Notes 3,4)
 - `dut1`   -- UT1-UTC (seconds, Note 5)
 - `elong`  -- longitude (radians, east +ve, Note 6)
 - `phi`    -- geodetic latitude (radians, Note 6)
 - `hm`     -- height above ellipsoid (m, geodetic Notes 6,8)
 - `xp, yp` -- polar motion coordinates (radians, Note 7)
 - `phpa`   -- pressure at the observer (hPa = mB, Note 8)
 - `tc`     -- ambient temperature at the observer (deg C)
 - `rh`     -- relative humidity at the observer (range 0-1)
 - `wl`     -- wavelength (micrometers, Note 9)

# Output

 - `rc, dc` -- ICRS astrometric RA,Dec (radians)

# Note

1) "Observed" Az,ZD means the position that would be seen by a perfect
   geodetically aligned theodolite.  (Zenith distance is used rather
   than altitude in order to reflect the fact that no allowance is
   made for depression of the horizon.)  This is related to the
   observed HA,Dec via the standard rotation, using the geodetic
   latitude (corrected for polar motion), while the observed HA and RA
   are related simply through the Earth rotation angle and the site
   longitude.  "Observed" RA,Dec or HA,Dec thus means the position
   that would be seen by a perfect equatorial with its polar axis
   aligned to the Earth's axis of rotation.

2) Only the first character of the type argument is significant.  "R"
   or "r" indicates that ob1 and ob2 are the observed right ascension
   and declination; "H" or "h" indicates that they are hour angle
   (west +ve) and declination; anything else ("A" or "a" is
   recommended) indicates that ob1 and ob2 are azimuth (north zero,
   east 90 deg) and zenith distance.

3) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
   convenient way between the two arguments, for example where utc1 is
   the Julian Day Number and utc2 is the fraction of a day.

   However, JD cannot unambiguously represent UTC during a leap second
   unless special measures are taken.  The convention in the present
   function is that the JD day represents UTC days whether the length
   is 86399, 86400 or 86401 SI seconds.

   Applications should use the function eraDtf2d to convert from
   calendar date and time of day into 2-part quasi Julian Date, as it
   implements the leap-second-ambiguity convention just described.

4) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale or that are too far in the future to
   be trusted.  See eraDat for further details.

5) UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
   one second at the end of each positive UTC leap second, introduced
   in order to keep UT1-UTC within +/- 0.9s.  n.b. This practice is
   under review, and in the future UT1-UTC may grow essentially
   without limit.

6) The geographical coordinates are with respect to the ERFA_WGS84
   reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN: the
   longitude required by the present function is east-positive
   (i.e. right-handed), in accordance with geographical convention.

7) The polar motion xp,yp can be obtained from IERS bulletins.  The
   values are the coordinates (in radians) of the Celestial
   Intermediate Pole with respect to the International Terrestrial
   Reference System (see IERS Conventions 2003), measured along the
   meridians 0 and 90 deg west respectively.  For many applications,
   xp and yp can be set to zero.

8) If hm, the height above the ellipsoid of the observing station in
   meters, is not known but phpa, the pressure in hPa (=mB), is
   available, an adequate estimate of hm can be obtained from the
   expression

		  hm = -29.3 * tsl * log ( phpa / 1013.25 );

   where tsl is the approximate sea-level air temperature in K (See
   Astrophysical Quantities, C.W.Allen, 3rd edition, section 52).
   Similarly, if the pressure phpa is not known, it can be estimated
   from the height of the observing station, hm, as follows:

		  phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

   Note, however, that the refraction is nearly proportional to the
   pressure and that an accurate phpa value is important for precise
   work.

9) The argument wl specifies the observing wavelength in micrometers.
   The transition from optical to radio is assumed to occur at 100
   micrometers (about 3000 GHz).

10) The accuracy of the result is limited by the corrections for
	refraction, which use a simple A*tan(z) + B*tan^3(z) model.
	Providing the meteorological parameters are known accurately and
	there are no gross local effects, the predicted astrometric
	coordinates should be within 0.05 arcsec (optical) or 1 arcsec
	(radio) for a zenith distance of less than 70 degrees, better than
	30 arcsec (optical or radio) at 85 degrees and better than 20
	arcmin (optical) or 30 arcmin (radio) at the horizon.

	Without refraction, the complementary functions eraAtco13 and
	eraAtoc13 are self-consistent to better than 1 microarcsecond all
	over the celestial sphere.  With refraction included, consistency
	falls off at high zenith distances, but is still better than 0.05
	arcsec at 85 degrees.

11) It is advisable to take great care with units, as even unlikely
	values of the input parameters are accepted and processed in
	accordance with the models used.
"""
function atoc13(tp::Char, ob1::AbstractFloat, ob2::AbstractFloat, utc1::AbstractFloat, utc2::AbstractFloat,
	dut1::AbstractFloat, elong::AbstractFloat, ϕ::AbstractFloat, hm::AbstractFloat, xp::AbstractFloat,
	yp::AbstractFloat, phpa::AbstractFloat, tc::AbstractFloat, rh::AbstractFloat, wl::AbstractFloat)
	#  Star-independent astrometry parameters
	a, eo = apco13(utc1, utc2, dut1, elong, ϕ, hm, xp, yp, phpa, tc, rh, wl)
	#  Transform observed to CIRS
	ri, di = atoiq(tp, ob1, ob2, a)
	#  Transform CIRS to ICRS
	ra, dec = values(aticq(ri, di, a))
	(ra = ra, dec = dec)
end

"""
	atoi13(tp::Char, ob1::AbstractFloat, ob2::AbstractFloat, utc1::AbstractFloat, utc2::AbstractFloat,
		   dut1::AbstractFloat, elong::AbstractFloat, ϕ::AbstractFloat, hm::AbstractFloat, xp::AbstractFloat,
		   yp::AbstractFloat, phpa::AbstractFloat, tc::AbstractFloat, rh::AbstractFloat, wl::AbstractFloat)

Observed place to CIRS.  The caller supplies UTC, site coordinates,
ambient air conditions and observing wavelength.

# Input

 - `type`   -- type of coordinates - "R", "H" or "A" (Notes 1,2)
 - `ob1`    -- observed Az, HA or RA (radians; Az is N=0,E=90)
 - `ob2`    -- observed ZD or Dec (radians)
 - `utc1`   -- UTC as a 2-part...
 - `utc2`   -- ...quasi Julian Date (Notes 3,4)
 - `dut1`   -- UT1-UTC (seconds, Note 5)
 - `elong`  -- longitude (radians, east +ve, Note 6)
 - `phi`    -- geodetic latitude (radians, Note 6)
 - `hm`     -- height above the ellipsoid (meters, Notes 6,8)
 - `xp, yp` -- polar motion coordinates (radians, Note 7)
 - `phpa`   -- pressure at the observer (hPa = mB, Note 8)
 - `tc`     -- ambient temperature at the observer (deg C)
 - `rh`     -- relative humidity at the observer (range 0-1)
 - `wl`     -- wavelength (micrometers, Note 9)

# Output

 - `ri`     -- CIRS right ascension (CIO-based, radians)
 - `di`     -- CIRS declination (radians)

# Note

1) "Observed" Az,ZD means the position that would be seen by a perfect
   geodetically aligned theodolite.  (Zenith distance is used rather
   than altitude in order to reflect the fact that no allowance is
   made for depression of the horizon.)  This is related to the
   observed HA,Dec via the standard rotation, using the geodetic
   latitude (corrected for polar motion), while the observed HA and RA
   are related simply through the Earth rotation angle and the site
   longitude.  "Observed" RA,Dec or HA,Dec thus means the position
   that would be seen by a perfect equatorial with its polar axis
   aligned to the Earth's axis of rotation.

2) Only the first character of the type argument is significant.  "R"
   or "r" indicates that ob1 and ob2 are the observed right ascension
   and declination; "H" or "h" indicates that they are hour angle
   (west +ve) and declination; anything else ("A" or "a" is
   recommended) indicates that ob1 and ob2 are azimuth (north zero,
   east 90 deg) and zenith distance.

3) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
   convenient way between the two arguments, for example where utc1 is
   the Julian Day Number and utc2 is the fraction of a day.

   However, JD cannot unambiguously represent UTC during a leap second
   unless special measures are taken.  The convention in the present
   function is that the JD day represents UTC days whether the length
   is 86399, 86400 or 86401 SI seconds.

   Applications should use the function eraDtf2d to convert from
   calendar date and time of day into 2-part quasi Julian Date, as it
   implements the leap-second-ambiguity convention just described.

4) The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale or that are too far in the future to
   be trusted.  See eraDat for further details.

5) UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
   one second at the end of each positive UTC leap second, introduced
   in order to keep UT1-UTC within +/- 0.9s.  n.b. This practice is
   under review, and in the future UT1-UTC may grow essentially
   without limit.

6) The geographical coordinates are with respect to the ERFA_WGS84
   reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN: the
   longitude required by the present function is east-positive
   (i.e. right-handed), in accordance with geographical convention.

7) The polar motion xp,yp can be obtained from IERS bulletins.  The
   values are the coordinates (in radians) of the Celestial
   Intermediate Pole with respect to the International Terrestrial
   Reference System (see IERS Conventions 2003), measured along the
   meridians 0 and 90 deg west respectively.  For many applications,
   xp and yp can be set to zero.

8) If hm, the height above the ellipsoid of the observing station in
   meters, is not known but phpa, the pressure in hPa (=mB), is
   available, an adequate estimate of hm can be obtained from the
   expression

		  hm = -29.3 * tsl * log ( phpa / 1013.25 );

   where tsl is the approximate sea-level air temperature in K (See
   Astrophysical Quantities, C.W.Allen, 3rd edition, section 52).
   Similarly, if the pressure phpa is not known, it can be estimated
   from the height of the observing station, hm, as follows:

		  phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );

   Note, however, that the refraction is nearly proportional to the
   pressure and that an accurate phpa value is important for precise
   work.

9) The argument wl specifies the observing wavelength in micrometers.
   The transition from optical to radio is assumed to occur at 100
   micrometers (about 3000 GHz).

10) The accuracy of the result is limited by the corrections for
	refraction, which use a simple A*tan(z) + B*tan^3(z) model.
	Providing the meteorological parameters are known accurately and
	there are no gross local effects, the predicted astrometric
	coordinates should be within 0.05 arcsec (optical) or 1 arcsec
	(radio) for a zenith distance of less than 70 degrees, better than
	30 arcsec (optical or radio) at 85 degrees and better than 20
	arcmin (optical) or 30 arcmin (radio) at the horizon.

	Without refraction, the complementary functions eraAtio13 and
	eraAtoi13 are self-consistent to better than 1 microarcsecond all
	over the celestial sphere.  With refraction included, consistency
	falls off at high zenith distances, but is still better than 0.05
	arcsec at 85 degrees.

12) It is advisable to take great care with units, as even unlikely
	values of the input parameters are accepted and processed in
	accordance with the models used.
"""
function atoi13(tp::Char, ob1::AbstractFloat, ob2::AbstractFloat, utc1::AbstractFloat, utc2::AbstractFloat,
	dut1::AbstractFloat, elong::AbstractFloat, ϕ::AbstractFloat, hm::AbstractFloat, xp::AbstractFloat,
	yp::AbstractFloat, phpa::AbstractFloat, tc::AbstractFloat, rh::AbstractFloat, wl::AbstractFloat)
	#  Star-independent astrometry parameters for CIRS->observed.
	#  Transform observed to CIRS.
	ri, di = atoiq(tp, ob1, ob2, apio13(utc1, utc2, dut1, elong, ϕ,
		hm, xp, yp, phpa, tc, rh, wl))
	(ra = ri, dec = di)
end

"""
	atoiq(tp::Char, ob1::AbstractFloat, ob2::AbstractFloat, a::Astrom)

Quick observed place to CIRS, given the star-independent astrometry
parameters.

Use of this function is appropriate when efficiency is important and
where many star positions are all to be transformed for one date.  The
star-independent astrometry parameters can be obtained by calling
eraApio[13] or eraApco[13].

# Input

 - `type`   -- type of coordinates: "R", "H" or "A" (Note 1)
 - `ob1`    -- observed Az, HA or RA (radians; Az is N=0,E=90)
 - `ob2`    -- observed ZD or Dec (radians)
 - `astrom` -- star-independent astrometry parameters:

# Output

 - `ri`     -- CIRS right ascension (CIO-based, radians)
 - `di`     -- CIRS declination (radians)

# Note

1) "Observed" Az,Zd means the position that would be seen by a perfect
   geodetically aligned theodolite.  This is related to the observed
   HA,Dec via the standard rotation, using the geodetic latitude
   (corrected for polar motion), while the observed HA and RA are
   related simply through the Earth rotation angle and the site
   longitude.  "Observed" RA,Dec or HA,Dec thus means the position
   that would be seen by a perfect equatorial with its polar axis
   aligned to the Earth's axis of rotation.  By removing from the
   observed place the effects of atmospheric refraction and diurnal
   aberration, the CIRS RA,Dec is obtained.

2) Only the first character of the type argument is significant.  "R"
   or "r" indicates that ob1 and ob2 are the observed right ascension
   and declination; "H" or "h" indicates that they are hour angle
   (west +ve) and declination; anything else ("A" or "a" is
   recommended) indicates that ob1 and ob2 are azimuth (north zero,
   east 90 deg) and zenith distance.  (Zenith distance is used rather
   than altitude in order to reflect the fact that no allowance is
   made for depression of the horizon.)

3) The accuracy of the result is limited by the corrections for
   refraction, which use a simple A*tan(z) + B*tan^3(z) model.
   Providing the meteorological parameters are known accurately and
   there are no gross local effects, the predicted intermediate
   coordinates should be within 0.05 arcsec (optical) or 1 arcsec
   (radio) for a zenith distance of less than 70 degrees, better than
   30 arcsec (optical or radio) at 85 degrees and better than 20
   arcmin (optical) or 25 arcmin (radio) at the horizon.

   Without refraction, the complementary functions eraAtioq and
   eraAtoiq are self-consistent to better than 1 microarcsecond all
   over the celestial sphere.  With refraction included, consistency
   falls off at high zenith distances, but is still better than 0.05
   arcsec at 85 degrees.

4) It is advisable to take great care with units, as even unlikely
   values of the input parameters are accepted and processed in
   accordance with the models used.
"""
function atoiq(tp::Char, ob1::AbstractFloat, ob2::AbstractFloat, a::Astrom)
	#  Minimum sin(alt) for refraction.
	SELMIN = 0.05
	#  if Az, ZD, convert to cartesian (S=0, E=90).
	aeo = MVector(0.0, 0.0, 0.0)
	if uppercase(tp) == 'A'
		aeo .= SVector(-cos(ob1)*sin(ob2), sin(ob1)*sin(ob2), cos(ob2))
	else
		#  If Ra, Dec, convert to HA, Dec.
		if uppercase(tp) == 'R'
			ob1 = a.eral - ob1
		end
		#  To cartesian -HA, Dec and then to cartesian Az, El (S=0, E=90).
		aeo .= SMatrix{3, 3}([a.sphi 0.0 -a.cphi; 0.0 1.0 0.0; a.cphi 0.0 a.sphi])*s2c(-ob1, ob2)
	end
	#  Azimuth (S=0, E=90).
	az = aeo[1] != 0.0 || aeo[2] != 0.0 ? atan(aeo[2], aeo[1]) : 0.0
	#  Sine of observed ZD, and observed ZD.
	zdo = atan(norm(aeo[1:2]), aeo[3])

	####    Refraction    ####
	#  Fast algorithm using two constant model.
	tz = norm(aeo[1:2])/maximum((aeo[3], SELMIN))
	zdt = zdo + a.refa*tz + a.refb*tz^3
	#  To cartesian AZ, ZD.
	aet = SVector(cos(az)*sin(zdt), sin(az)*sin(zdt), cos(zdt))
	#  Cartesian Az, ZD to cartesian -HA, Dec.
	mhda = SMatrix{3, 3}([a.sphi 0.0 a.cphi; 0.0 1.0 0.0; -a.cphi 0.0 a.sphi])*aet
	#  Diurnal aberration.
	hd = (1.0 + a.diurab*mhda[2]) .* (mhda .- SVector(0.0, a.diurab, 0.0))
	#  Polar motion.
	hma, dec = c2s(SMatrix{3, 3}(
		[cos(a.xpl)  sin(a.xpl)*sin(a.ypl) -sin(a.xpl)*cos(a.ypl);
			   0.0              cos(a.ypl)             sin(a.ypl);
			sin(a.xpl) -cos(a.xpl)*sin(a.ypl)  cos(a.xpl)*cos(a.ypl)])*hd)
	(ra = mod2pi(a.eral + hma), dec = dec)
end

"""
	ld(bm::AbstractFloat, p::AbstractVector{<:AbstractFloat}, q::AbstractVector{<:AbstractFloat}, e::AbstractVector{<:AbstractFloat},
	   em::AbstractFloat, dlim::AbstractFloat)

Apply light deflection by a solar-system body, as part of transforming
coordinate direction into natural direction.

# Input

 - `bm`     -- mass of the gravitating body (solar masses)
 - `p`      -- direction from observer to source (unit vector)
 - `q`      -- direction from body to source (unit vector)
 - `e`      -- direction from body to observer (unit vector)
 - `em`     -- distance from body to observer (AU)
 - `dlim`   -- deflection limiter (Note 4)

# Output

 - `p1`     -- observer to deflected source (unit vector)

# Note

1) The algorithm is based on Expr. (70) in Klioner (2003) and
   Expr. (7.63) in the Explanatory Supplement (Urban & Seidelmann
   2013), with some rearrangement to minimize the effects of machine
   precision.

2) The mass parameter bm can, as required, be adjusted in order to
   allow for such effects as quadrupole field.

3) The barycentric position of the deflecting body should ideally
   correspond to the time of closest approach of the light ray to the
   body.

4) The deflection limiter parameter dlim is phi^2/2, where phi is the
   angular separation (in radians) between source and body at which
   limiting is applied.  As phi shrinks below the chosen threshold,
   the deflection is artificially reduced, reaching zero for phi = 0.

5) The returned vector p1 is not normalized, but the consequential
   departure from unit magnitude is always negligible.

6) The arguments p and p1 can be the same array.

7) To accumulate total light deflection taking into account the
   contributions from several bodies, call the present function for
   each body in succession, in decreasing order of distance from the
   observer.

8) For efficiency, validation is omitted.  The supplied vectors must
   be of unit magnitude, and the deflection limiter non-zero and
   positive.

# References

Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to the
Astronomical Almanac, 3rd ed., University Science Books (2013).

Klioner, Sergei A., "A practical relativistic model for micro-
arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).
"""
function ld(bm::F, p::V, q::W, e::X, em::F, dlim::F) where
{F <: AbstractFloat, V <: AbstractVector{<:AbstractFloat},
	W <: AbstractVector{<:AbstractFloat}, X <: AbstractVector{<:AbstractFloat}}
	#  2*G*bm/(em*c^2*(q*(q+e))).
	#  Apply the deflection.
	p .+ SCHWARZRADIUS*bm/em/maximum((dot(q, (q .+ e)), dlim)) .* pxp(p, pxp(e, q))
end

"""
	ldn(n::Int, b::AbstractVector{Ldbody}, ob::AbstractVector{<:AbstractFloat}, sc::AbstractVector{<:AbstractFloat})

For a star, apply light deflection by multiple solar-system bodies, as
part of transforming coordinate direction into natural direction.

# Input

 - `n`     -- number of bodies (note 1)
 - `b`     -- data for each of the n bodies (Notes 1,2):
 - `bm`    -- mass of the body (solar masses, Note 3)
 - `dl`    -- deflection limiter (Note 4)
 - `pv`    -- barycentric PV of the body (AU, AU/day)
 - `ob`    -- barycentric position of the observer (AU)
 - `sc`    -- observer to star coord direction (unit vector)

# Output

 - `sn`    -- observer to deflected star (unit vector)

# Note

1) The array b contains n entries, one for each body to be considered.
   If n = 0, no gravitational light deflection will be applied, not
   even for the Sun.

2) The array b should include an entry for the Sun as well as for any
   planet or other body to be taken into account.  The entries should
   be in the order in which the light passes the body.

3) In the entry in the b array for body i, the mass parameter b[i].bm
   can, as required, be adjusted in order to allow for such effects as
   quadrupole field.

4) The deflection limiter parameter b[i].dl is phi^2/2, where phi is
   the angular separation (in radians) between star and body at which
   limiting is applied.  As phi shrinks below the chosen threshold,
   the deflection is artificially reduced, reaching zero for phi = 0.
   Example values suitable for a terrestrial observer, together with
   masses, are as follows:

	  body i     b[i].bm        b[i].dl

	  Sun        1.0            6e-6
	  Jupiter    0.00095435     3e-9
	  Saturn     0.00028574     3e-10

5) For cases where the starlight passes the body before reaching the
   observer, the body is placed back along its barycentric track by
   the light time from that point to the observer.  For cases where
   the body is "behind" the observer no such shift is applied.  If a
   different treatment is preferred, the user has the option of
   instead using the eraLd function.  Similarly, eraLd can be used for
   cases where the source is nearby, not a star.

6) The returned vector sn is not normalized, but the consequential
   departure from unit magnitude is always negligible.

7) The arguments sc and sn can be the same array.

8) For efficiency, validation is omitted.  The supplied masses must be
   greater than zero, the position and velocity vectors must be right,
   and the deflection limiter greater than zero.

# References

Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to the
Astronomical Almanac, 3rd ed., University Science Books (2013),
Section 7.2.4.
"""
function ldn(n::Int, b::AbstractVector{Ldbody}, ob::V, sc::W) where
{V <: AbstractVector{<:AbstractFloat}, W <: AbstractVector{<:AbstractFloat}}
	sn = sc[:]
	for body in b
		#  Body to observer vector at epoch of observation (AU).
		v = ob .- body.pv[1]
		#  Minus the time since the light passed the body (days).
		Δt = minimum((ASTRUNIT/LIGHTSPEED/SECPERDAY*dot(sn, v), 0.0))
		#  Backtrack the body to the time the light was passing the body.
		#  Body to observer vector as magnitude and direction.
		em, e = pn(v .- Δt .* body.pv[2])
		#  Apply light deflection for this body.
		sn .= ld(body.bm, sn, sn, e, em, body.dl)
	end
	sn
end

"""
	ldsun(p::AbstractVector{<:AbstractFloat}, e::AbstractVector{<:AbstractFloat}, em::AbstractFloat)

Deflection of starlight by the Sun.

# Input

 - `p`      -- direction from observer to star (unit vector)
 - `e`      -- direction from Sun to observer (unit vector)
 - `em`     -- distance from Sun to observer (AU)

# Output

 - `p1`     -- observer to deflected star (unit vector)

# Note

1) The source is presumed to be sufficiently distant that its
   directions seen from the Sun and the observer are essentially the
   same.

2) The deflection is restrained when the angle between the star and
   the center of the Sun is less than a threshold value, falling to
   zero deflection for zero separation.  The chosen threshold value is
   within the solar limb for all solar-system applications, and is
   about 5 arcminutes for the case of a terrestrial observer.

3) The arguments p and p1 can be the same array.
"""
function ldsun(p::V, e::W, em::F) where
{F <: AbstractFloat, V <: AbstractVector{<:AbstractFloat},
	W <: AbstractVector{<:AbstractFloat}}
	#  Deflection limiter (smaller for distant observers).
	#  Apply the deflection.
	ld(1.0, p, p, e, em, 1e-6/(em^2 > 1.0 ? em^2 : 1.0))
end

"""
	pmpx(rc::AbstractFloat, dc::AbstractFloat, pr::AbstractFloat, pd::AbstractFloat, px::AbstractFloat,
		 rv::AbstractFloat, pmt::AbstractFloat, pob::AbstractVector{<:AbstractFloat})

Proper motion and parallax.

# Input

 - `rc, dc` -- ICRS RA,Dec at catalog epoch (radians)
 - `pr`     -- RA proper motion (radians/year, Note 1)
 - `pd`     -- Dec proper motion (radians/year)
 - `px`     -- parallax (arcsec)
 - `rv`     -- radial velocity (km/s, +ve if receding)
 - `pmt`    -- proper motion time interval (SSB, Julian years)
 - `pob`    -- SSB to observer vector (AU)

# Output

 - `pco`    -- coordinate direction (BCRS unit vector)

# Note

1) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

2) The proper motion time interval is for when the starlight reaches
   the solar system barycenter.

3) To avoid the need for iteration, the Roemer effect (i.e. the small
   annual modulation of the proper motion coming from the changing
   light time) is applied approximately, using the direction of the
   star at the catalog epoch.

# References

1984 Astronomical Almanac, pp B39-B41.

Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to the
Astronomical Almanac, 3rd ed., University Science Books (2013),
Section 7.2.
"""
function pmpx(rc::AbstractFloat, dc::AbstractFloat, pr::AbstractFloat, pd::AbstractFloat,
	px::AbstractFloat, rv::AbstractFloat, pmt::AbstractFloat,
   pob::AbstractVector{<:AbstractFloat})
	#  Spherical coordinates to unit vector (and useful functions.)
	p = MVector(cos(rc)*cos(dc), sin(rc)*cos(dc), sin(dc))
	#  Space motion (radian per year).
	rvpx = SECPERDAY*(1000*DAYPERYEAR)/ASTRUNIT*rv*deg2rad(px/3600.0)
	pm = SVector(rvpx*p[1] - pr*p[2] - pd*cos(rc)*p[3],
		rvpx*p[2] + pr*p[1] - pd*sin(rc)*p[3],
		rvpx*p[3] + pd*cos(dc))
	#  Proper motion time interval (y) including Roemer effect.
	#  Coordinate direction of star (unit vector, BCRS).
	p .+= (pmt + AULIGHT*dot(p, pob)) .* pm - deg2rad(px/3600.0) .* pob
	p ./ norm(p)
end

"""
	pmsafe(ra::AbstractFloat, dec::AbstractFloat, pmr::AbstractFloat, pmd::AbstractFloat, px::AbstractFloat,
		   rv::AbstractFloat, ep1a::AbstractFloat, ep1b::AbstractFloat, ep2a::AbstractFloat,
		   ep2b::AbstractFloat)

Star proper motion: update star catalog data for space motion, with
special handling to handle the zero parallax case.

# Input

 - `ra1`    -- right ascension (radians), before
 - `dec1`   -- declination (radians), before
 - `pmr1`   -- RA proper motion (radians/year), before
 - `pmd1`   -- Dec proper motion (radians/year), before
 - `px1`    -- parallax (arcseconds), before
 - `rv1`    -- radial velocity (km/s, +ve = receding), before
 - `ep1a`   -- "before" epoch, part A (Note 1)
 - `ep1b`   -- "before" epoch, part B (Note 1)
 - `ep2a`   -- "after" epoch, part A (Note 1)
 - `ep2b`   -- "after" epoch, part B (Note 1)

# Output

 - `ra2`    -- right ascension (radians), after
 - `dec2`   -- declination (radians), after
 - `pmr2`   -- RA proper motion (radians/year), after
 - `pmd2`   -- Dec proper motion (radians/year), after
 - `px2`    -- parallax (arcseconds), after
 - `rv2`    -- radial velocity (km/s, +ve = receding), after

# Note

1) The starting and ending TDB epochs ep1a+ep1b and ep2a+ep2b are
   Julian Dates, apportioned in any convenient way between the two
   parts (A and B).  For example, JD(TDB)=2450123.7 could be expressed
   in any of these ways, among others:

		  epNa            epNb

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

2) In accordance with normal star-catalog conventions, the object's
   right ascension and declination are freed from the effects of
   secular aberration.  The frame, which is aligned to the catalog
   equator and equinox, is Lorentzian and centered on the SSB.

   The proper motions are the rate of change of the right ascension
   and declination at the catalog epoch and are in radians per TDB
   Julian year.

   The parallax and radial velocity are in the same frame.

3) Care is needed with units.  The star coordinates are in radians and
   the proper motions in radians per Julian year, but the parallax is
   in arcseconds.

4) The RA proper motion is in terms of coordinate angle, not true
   angle.  If the catalog uses arcseconds for both RA and Dec proper
   motions, the RA proper motion will need to be divided by cos(Dec)
   before use.

5) Straight-line motion at constant speed, in the inertial frame, is
   assumed.

6) An extremely small (or zero or negative) parallax is overridden to
   ensure that the object is at a finite but very large distance, but
   not so large that the proper motion is equivalent to a large but
   safe speed (about 0.1c using the chosen constant).  A warning
   status of 1 is added to the status if this action has been taken.

7) If the space velocity is a significant fraction of c (see the
   constant VMAX in the function eraStarpv), it is arbitrarily set to
   zero.  When this action occurs, 2 is added to the status.

8) The relativistic adjustment carried out in the eraStarpv function
   involves an iterative calculation.  If the process fails to
   converge within a set number of iterations, 4 is added to the
   status.
"""
function pmsafe(ra::AbstractFloat, dec::AbstractFloat, pmr::AbstractFloat, pmd::AbstractFloat,
	plx::AbstractFloat, rv::AbstractFloat, ep1a::AbstractFloat, ep1b::AbstractFloat,
	ep2a::AbstractFloat, ep2b::AbstractFloat)
	### println("$ra, $dec, $pmr, $pmd, $px, $rv, $ep1a, $ep1b, $ep2a, $ep2b")
	#  Minimum allowed parallax and factor giving maximum allowed transverse
	# speed of about 1% c.
	PXMIN, F = 5e-7, 326.0
	#  Proper motion in one year (radians)
	pm = F*seps(ra, dec, ra+pmr, dec+pmd)
	#  Override the parallax to reduce chances of a warning status.
	if plx < pm
		plx = pm      # @warn
	elseif plx < PXMIN
		plx = PXMIN   # @warn
	end
	#  Carry out the transformation using the modified parallax.
	### println("$ra, $dec, $pmr, $pmd, $px, $rv, $ep1a, $ep1b, $ep2a, $ep2b")
	starpm(ra, dec, pmr, pmd, plx, rv, ep1a, ep1b, ep2a, ep2b)
end

"""
	pvtob(elong::AbstractFloat, ϕ::AbstractFloat, hm::AbstractFloat, xp::AbstractFloat, yp::AbstractFloat,
		  sp::AbstractFloat, θ::AbstractFloat)

Position and velocity of a terrestrial observing station.

# Input

 - `elong`  -- longitude (radians, east +ve, Note 1)
 - `phi`    -- latitude (geodetic, radians, Note 1)
 - `hm`     -- height above ref. ellipsoid (geodetic, m)
 - `xp, yp` -- coordinates of the pole (radians, Note 2)
 - `sp`     -- the TIO locator s' (radians, Note 2)
 - `theta`  -- Earth rotation angle (radians, Note 3)

# Output

 - `pv`     -- position/velocity vector (m, m/s, CIRS)

# Note

1) The terrestrial coordinates are with respect to the ERFA_WGS84
   reference ellipsoid.

2) xp and yp are the coordinates (in radians) of the Celestial
   Intermediate Pole with respect to the International Terrestrial
   Reference System (see IERS Conventions), measured along the
   meridians 0 and 90 deg west respectively.  sp is the TIO locator
   s', in radians, which positions the Terrestrial Intermediate Origin
   on the equator.  For many applications, xp, yp and (especially) sp
   can be set to zero.

3) If theta is Greenwich apparent sidereal time instead of Earth
   rotation angle, the result is with respect to the true equator and
   equinox of date, i.e. with the x-axis at the equinox rather than
   the celestial intermediate origin.

4) The velocity units are meters per UT1 second, not per SI second.
   This is unlikely to have any practical consequences in the modern
   era.

5) No validation is performed on the arguments.  Error cases that
   could lead to arithmetic exceptions are trapped by the eraGd2gc
   function, and the result set to zeros.

# References

McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003), IERS
Technical Note No. 32, BKG (2004)

Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to the
Astronomical Almanac, 3rd ed., University Science Books (2013),
Section 7.4.3.3.
"""
function pvtob(elong::AbstractFloat, ϕ::AbstractFloat, hm::AbstractFloat, xp::AbstractFloat,
	yp::AbstractFloat, sp::AbstractFloat, θ::AbstractFloat)
	#  Earth rotation rate in radians per UT1 seconds.
	Ω = 2π*Ω_Earth_2003/SECPERDAY
	#  Geodetic to geocentric transformation (WGS84), polar motion and
	#  TIO position.
	x, y, z = trxp(pom00(xp, yp, sp), gd2gc(:WGS84, elong, ϕ, hm))
	#  Functions of ERA, position and velocity.
	pv = SVector(SVector(sum(sincos(θ) .* (-y, x)), sum(sincos(θ) .* (x, y)), z),
		Ω .* SVector(sum(sincos(θ) .* (-x, -y)), sum(sincos(θ) .* (-y, x)), 0.0))
end

"""
	refco(phpa::AbstractFloat, tc::AbstractFloat, rh::AbstractFloat, wl::AbstractFloat)

Determine the constants A and B in the atmospheric refraction model dZ
= A tan Z + B tan^3 Z.

Z is the "observed" zenith distance (i.e. affected by refraction) and
dZ is what to add to Z to give the "topocentric" (i.e. in vacuo)
zenith distance.

# Input

 - `phpa` -- pressure at the observer (hPa = millibar)
 - `tc`   -- ambient temperature at the observer (deg C)
 - `rh`   -- humidity at the observer (range 0-1)
 - `wl`   -- double    wavelength (micrometers)

# Output

 - `refa` -- tan Z coefficient (radians)
 - `refb` -- tan^3 Z coefficient (radians)

# Note

1) The model balances speed and accuracy to give good results in
   applications where performance at low altitudes is not paramount.
   Performance is maintained across a range of conditions, and applies
   to both optical/IR and radio.

2) The model omits the effects of (i) height above sea level (apart
   from the reduced pressure itself), (ii) latitude (i.e. the
   flattening of the Earth), (iii) variations in tropospheric lapse
   rate and (iv) dispersive effects in the radio.  The model was
   tested using the following range of conditions:

	 lapse rates 0.0055, 0.0065, 0.0075 deg/meter
	 latitudes 0, 25, 50, 75 degrees
	 heights 0, 2500, 5000 meters ASL
	 pressures mean for height -10% to +5% in steps of 5%
	 temperatures -10 deg to +20 deg with respect to 280 deg at SL
	 relative humidity 0, 0.5, 1
	 wavelengths 0.4, 0.6, ... 2 micron, + radio
	 zenith distances 15, 45, 75 degrees

   The accuracy with respect to raytracing through a model atmosphere
   was as follows:

						  worst         RMS

	 optical/IR           62 mas       8 mas
	 radio               319 mas      49 mas

   For this particular set of conditions:

	 lapse rate 0.0065 K/meter
	 latitude 50 degrees
	 sea level
	 pressure 1005 mb
	 temperature 280.15 K
	 humidity 80%
	 wavelength 5740 Angstroms

   the results were as follows:
	 ZD       raytrace     eraRefco   Saastamoinen

	 10         10.27        10.27        10.27
	 20         21.19        21.20        21.19
	 30         33.61        33.61        33.60
	 40         48.82        48.83        48.81
	 45         58.16        58.18        58.16
	 50         69.28        69.30        69.27
	 55         82.97        82.99        82.95
	 60        100.51       100.54       100.50
	 65        124.23       124.26       124.20
	 70        158.63       158.68       158.61
	 72        177.32       177.37       177.31
	 74        200.35       200.38       200.32
	 76        229.45       229.43       229.42
	 78        267.44       267.29       267.41
	 80        319.13       318.55       319.10

	deg        arcsec       arcsec       arcsec

   The values for Saastamoinen's formula (which includes terms up to
   tan^5) are taken from Hohenkerk and Sinclair (1985).

3) A wl value in the range 0-100 selects the optical/IR case and is
   wavelength in micrometers.  Any value outside this range selects
   the radio case.

4) Outlandish input parameters are silently limited to mathematically
   safe values.  Zero pressure is permissible, and causes zeroes to be
   returned.

5) The algorithm draws on several sources, as follows:

   a) The formula for the saturation vapour pressure of water as a
	  function of temperature and temperature is taken from Equations
	  (A4.5-A4.7) of Gill (1982).

   b) The formula for the water vapour pressure, given the saturation
	  pressure and the relative humidity, is from Crane (1976),
	  Equation (2.5.5).

   c) The refractivity of air is a function of temperature, total
	  pressure, water-vapour pressure and, in the case of optical/IR,
	  wavelength.  The formulae for the two cases are developed from
	  Hohenkerk & Sinclair (1985) and Rueger (2002).  The IAG (1999)
	  optical refractivity for dry air is used.

   d) The formula for beta, the ratio of the scale height of the
	  atmosphere to the geocentric distance of the observer, is an
	  adaption of Equation (9) from Stone (1996).  The adaptations,
	  arrived at empirically, consist of (i) a small adjustment to the
	  coefficient and (ii) a humidity term for the radio case only.

   e) The formulae for the refraction constants as a function of n-1
	  and beta are from Green (1987), Equation (4.31).

# References

Crane, R.K., Meeks, M.L. (ed), "Refraction Effects in the Neutral
Atmosphere", Methods of Experimental Physics: Astrophysics 12B,
Academic Press, 1976.

Gill, Adrian E., "Atmosphere-Ocean Dynamics", Academic Press, 1982.

Green, R.M., "Spherical Astronomy", Cambridge University Press, 1987.

Hohenkerk, C.Y., & Sinclair, A.T., NAO Technical Note No. 63, 1985.

IAG Resolutions adopted at the XXIIth General Assembly in Birmingham,
1999, Resolution 3.

Rueger, J.M., "Refractive Index Formulae for Electronic Distance
Measurement with Radio and Millimetre Waves", in Unisurv Report S-68,
School of Surveying and Spatial Information Systems, University of New
South Wales, Sydney, Australia, 2002.

Stone, Ronald C., P.A.S.P. 108, 1051-1058, 1996.
"""
function refco(phpa::AbstractFloat, tc::AbstractFloat, rh::AbstractFloat, wl::AbstractFloat)
	#  Determine the spectral band, either optical/IR or radio, based on
	#  the wavelength of > or < 100 microns.
	optical = wl <= 100.0 ? true : false
	#  Restrict parameters to safe values.
	t = minimum((maximum((tc, -150.0)), 200.0))
	p = minimum((maximum((phpa, 0.0)), 10000.0))
	r = minimum((maximum((rh, 0.0)), 1.0))
	w = minimum((maximum((wl, 0.1)), 1e6))
	#  Water vapor pressure at the observer.
	pw = p > 0.0 ? h20pres(p, t, r) : 0.0
	#  Refractive index minus 1 at the observer.
	if optical
		wpco = (77.53484e-6, 4.39108e-7, 3.666e-9)
		γ = (Polynomial(wpco...)(1/w^2)*p - 11.2684e-6*pw)/(273.15 + t)
	else
		γ = (77.6890e-6*p - (6.3938e-6 - 0.375463/(273.15 + t))*pw)/(273.15 + t)
	end
	#  Formula for β from Stone, with empirical adjustments.
	β = 4.4474e-6*(273.15 + t)
	β -= !optical ? 0.0074*pw*β : 0.0
	#  Refraction constants from Green.
	refa, refb = γ*(1.0 - β), -γ*(β - γ/2.0)
end

"""
	h20pres(p::AbstractFloat, t::AbstractFloat, r::AbstractFloat)

Calculate the water vapor pressure given temperature and pressure.

# Input

 - `p`  -- atmospheric pressure
 - `t`  -- atmospheric temperature
 - `r`  -- relative humidity

# Output

 - `wp` -- water vapor pressure

"""
function h20pres(p::AbstractFloat, t::AbstractFloat, r::AbstractFloat)
	ps = (1.0 + p*(4.5e-6 + 6e-10*t^2)) *
		 10.0^((0.7859 + 0.03477*t)/(1.0 + 0.00412*t))
	r*ps/(1.0 - (1.0 - r)*ps/p)
end
