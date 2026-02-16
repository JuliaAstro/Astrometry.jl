#### Vector - Matrix / Angle Operations

"""
    a2af(ndp::Integer, angle::Real)

Decompose radians into degrees, arcminutes, arcseconds, and fraction.

# Input

 - `npd`   --number of useful digits
 - `angle` -- angle in radians

# Output

 - `dms`   -- angle in sign, degrees, minutes, seconds, and fraction

# Note

1) The argument ndp is interpreted as follows:

   ndp         resolution
    :      ...0000 00 00
   -7         1000 00 00
   -6          100 00 00
   -5           10 00 00
   -4            1 00 00
   -3            0 10 00
   -2            0 01 00
   -1            0 00 10
    0            0 00 01
    1            0 00 00.1
    2            0 00 00.01
    3            0 00 00.001
    :            0 00 00.000...

2) The largest positive useful value for ndp is determined by the size
   of angle, the format of Float64 on the target platform, and the
   risk of overflowing dms[3].  On a typical platform, for angle up to
   2pi, the available floating-point precision might correspond to
   ndp=12.  However, the practical limit is typically ndp=9, set by
   the capacity of a 32-bit int, or ndp=4 if int is only 16 bits.

3) The absolute value of angle may exceed 2pi.  In cases where it does
   not, it is up to the caller to test for and handle the case where
   days is very nearly 2pi and rounds up to 2pi and rounds up to 360
   degree, by testing for dms[0]=360 and setting dms[0-3] to zero.
"""
function a2af(ndp::Integer, angle::Real)
    NamedTuple{(:sign, :degree, :minute, :second, :fraction)}
    (values(d2tf(ndp, angle*15/2/pi)))
end

"""
    a2tf(ndp::Integer, angle::Real)

Decompose radians into hours, minutes, seconds, and fraction.

# Input

 - `npd`   -- number of useful digits
 - `angle` -- angle in radians

# Output
 - `hms`   -- angle in sign, hour, minutes, seconds, and fraction

# Note

1) The argument ndp is interpreted as follows:

   ndp         resolution
    :      ...0000 00 00
   -7         1000 00 00
   -6          100 00 00
   -5           10 00 00
   -4            1 00 00
   -3            0 10 00
   -2            0 01 00
   -1            0 00 10
    0            0 00 01
    1            0 00 00.1
    2            0 00 00.01
    3            0 00 00.001
    :            0 00 00.000...

2) The largest positive useful value for ndp is determined by the size
   of angle, the format of Float64 on the target platform, and the
   risk of overflowing hms[3].  On a typical platform, for angle up to
   2pi, the available floating-point precision might correspond to
   ndp=12.  However, the practical limit is typically ndp=9, set by
   the capacity of a 32-bit int, or ndp=4 if int is only 16 bits.

3) The absolute value of angle may exceed 2pi.  In cases where it does
   not, it is up to the caller to test for and handle the case where
   days is very nearly 2pi and rounds up to 2pi and rounds up to 24
   hours, by testing for ihmsf[0]=24 and setting ihmsf[0-3] to zero.
"""
function a2tf(ndp::Integer, angle::Real)
    NamedTuple{(:sign, :hour, :minute, :second, :fraction)}(d2tf(ndp, angle/2/pi))
end

"""
    af2a(sign::Char, degree::Integer, minute::Integer, second::Real)

Convert degrees, arcminutes, arcseconds to radians.

# Input

 - `sign`   -- sign of arc
 - `degree` -- degrees of arc
 - `minute` -- minutes of arc
 - `second` -- seconds of arc

# Output

 - `angle`  -- angle in radians

# Note

1) The result is computed even if any of the range checks fail.

2) Negative ideg, iamin and/or asec produce a warning status, but the
   absolute value is used in the conversion.

3) If there are multiple errors, the status value reflects only the
   first, the smallest taking precedence.
"""
function af2a(sign::Char, degree::Integer, minute::Integer, second::Real)
    @assert 0   <= degree < 360   "degree out of range [0-359]."
    @assert 0   <= minute <  60   "minute out of range [0-59]."
    @assert 0.0 <= second <  60.0 "second out of range [0-60]."
    
    onet = one(typeof(second))
    deg2rad(1/3600) * (sign == '-' ? -onet : onet) *
       (60*onet * (60*onet * abs(degree) + abs(minute)) + abs(second))
end

"""
    anp(angle::Real)

Normalize angle into the range 0 <= a < 2p.

# Input

 - `angle` -- angle in radians

# Output

 - `angle` -- angle in radians in range 0-2pi
"""
function anp(angle::Real)
    mod2pi(angle)
end

"""
    anpm(angle::Real)

Normalize angle into the range -pi <= a < +pi

# Input

 - `angle` -- angle in radians

# Output

 - `angle` -- angle in radians in range +/-pi
"""
function anpm(angle::Real)
    rem2pi(angle, RoundNearest)
end

"""
    d2tf(ndp::Integer, day::Real)

Decompose days to sign, hours, minutes, seconds, fraction.

# Input

 - `npd`   -- number of usefule digits
 - `day`   -- interval in days

# Output

 - `hms`   -- hms in sign, hours, minutes, seconds, and fraction

# Note

1) The argument ndp is interpreted as follows:

   ndp         resolution
    :      ...0000 00 00
   -7         1000 00 00
   -6          100 00 00
   -5           10 00 00
   -4            1 00 00
   -3            0 10 00
   -2            0 01 00
   -1            0 00 10
    0            0 00 01
    1            0 00 00.1
    2            0 00 00.01
    3            0 00 00.001
    :            0 00 00.000...

2) The largest positive useful value for ndp is determined by the size
   of days, the format of Float64 on the target platform, and the risk
   of overflowing hms[3].  On a typical platform, for days up to 1.0,
   the available floating-point precision might correspond to ndp=12.
   However, the practical limit is typically ndp=9, set by the
   capacity of a 32-bit int, or ndp=4 if int is only 16 bits.

3) The absolute value of days may exceed 1.0.  In cases where it does
   not, it is up to the caller to test for and handle the case where
   days is very nearly 1.0 and rounds up to 24 hours, by testing for
   ihmsf[0]=24 and setting ihmsf[0-3] to zero.
"""
function d2tf(ndp::Integer, day::Real)

    a = SECPERDAY * abs(day)
    if ndp < 0
        rs = prod(n == 2 || n == 4 ? 6 : 10 for n=1:-ndp)
        a = rs * round(a/rs)
    else
        rs = 10^ndp
    end

    rh, rm = 3600.0 * rs, 60.0 * rs

    sn = day >= 0.0 ? '+' : '-'
    a  = round(rs * a)
    ah = convert(Integer, trunc(a/rh))
    a -= ah*rh
    am = convert(Integer, trunc(a/rm))
    a -= am*rm
    as = convert(Integer, trunc(a/rs))
    af = convert(Integer, a - as*rs)

    NamedTuple{(:sign, :hour, :minute, :second, :fraction)}((sn, ah, am, as, af))
end

"""
    tf2a(sign::Char, hour::Integer, minute::Integer, second::Real)

Convert hours, minutes, seconds to radians.

# Input

 - `sign`   -- sign:  '-' = negative, otherwise positive
 - `hour`   -- hours
 - `minute` -- minutes
 - `second` -- seconds

# Output

 - `angle`  -- angle in radians

# Note

1) The result is computed even if any of the range checks fail.

2) Negative ihour, imin and/or sec produce a warning status, but the
   absolute value is used in the conversion.

3) If there are multiple errors, the status value reflects only the
   first, the smallest taking precedence.
"""
function tf2a(sign::Char, hour::Integer, minute::Integer, second::Real)
    @assert 0   <= hour   <   24  "hour out of range [0-23]."
    @assert 0   <= minute <   60  "minute out of range [0-59]."
    @assert 0.0 <= second < 60.0 "second out of range [0-60]."
    
    15*deg2rad(1/3600) * (sign == '-' ? -1.0 : 1.0) *
        (60.0 * (60.0 * abs(hour) + abs(minute)) + abs(second))
end

"""
    tf2d(sign::Char, hour::Integer, minute::Integer, second::Real)

Convert hours, minutes, seconds to days.

# Input

 - `sign`   -- sign:  '-' = negative, otherwise positive
 - `hour`   -- hours
 - `minute` -- minutes
 - `second` -- seconds

# Output

 - `day`    -- interval in days

# Note

1) The result is computed even if any of the range checks fail.

2) Negative ihour, imin and/or sec produce a warning status, but the
   absolute value is used in the conversion.

3) If there are multiple errors, the status value reflects only the
   first, the smallest taking precedence.
"""
function tf2d(sign::Char, hour::Integer, minute::Integer, second::Real)
    @assert 0   <= hour   <   24  "hour out of range [0-23]."
    @assert 0   <= minute <   60  "minute out of range [0-59]."
    @assert 0.0 <= second < 60.0 "second out of range [0-60]."
    
    (sign == '-' ? -1.0 : 1.0) *
        (60.0 * (60.0 * abs(hour) + abs(minute)) + abs(second)) / SECPERDAY
end

#### Vector - Matrix / Build Rotations

"""
    rx(ϕ::Real, r::AbstractMatrix{<:Real})

Rotate an r-matrix about the x-axis.

# Input

 - `ϕ`     -- angle (radians)
 - `r`     -- r-matrix

# Output

 - `r`     -- r-matrix, rotated

# Note

1) Calling this function with positive ϕ incorporates in the supplied
   r-matrix r an additional rotation, about the x-axis, anticlockwise
   as seen looking towards the origin from positive x.

2) The additional rotation can be represented by this matrix:

       (  1        0            0      )
       (                               )
       (  0   + cos(ϕ)   + sin(ϕ)  )
       (                               )
       (  0   - sin(ϕ)   + cos(ϕ)  )

"""
function rx(ϕ::Real, r::AbstractMatrix{<:Real})
    #  Matrix multiplication performs two allocations
    Rx(ϕ)*r
end

"""
    ry(θ::Real, r::AbstractMatrix{<:Real})

Rotate an r-matrix about the y-axis.

# Input

 - `θ`     -- angle (radians)
 - `r`     -- r-matrix

# Output

 - `r`     -- r-matrix, rotated

# Note

1) Calling this function with positive theta incorporates in the
   supplied r-matrix r an additional rotation, about the y-axis,
   anticlockwise as seen looking towards the origin from positive y.

2) The additional rotation can be represented by this matrix:

       (  + cos(θ)     0      - sin(θ)  )
       (                                        )
       (       0           1           0        )
       (                                        )
       (  + sin(θ)     0      + cos(θ)  )
"""
function ry(θ::Real, r::AbstractMatrix{<:Real})
    #  Matrix multiplication performs two allocations
    Ry(θ)*r
end

"""
    rz(ψ::Real, r::AbstractMatrix{<:Real})

Rotate an r-matrix about the z-axis.

# Input

 - `ψ`     -- angle (radians)
 - `r `    -- r-matrix

# Output

 - `r`     -- r-matrix, rotated

# Note

1) Calling this function with positive ψ incorporates in the supplied
   r-matrix r an additional rotation, about the z-axis, anticlockwise
   as seen looking towards the origin from positive z.

2) The additional rotation can be represented by this matrix:

       (  + cos(ψ)   + sin(ψ)     0  )
       (                              )
       (  - sin(ψ)   + cos(ψ)     0  )
       (                              )
       (       0            0      1  )

"""
function rz(ψ::Real, r::AbstractMatrix{<:Real})
    #  Matrix multiplication performs two allocations
    Rz(ψ)*r
end

#### Vector - Matrix / Copy, Extend, Extract

"""
    cp(p::AbstractVector{<:Real})

Copy a p-vector.

# Input

 - `p`     -- p-vector to be copied

# Output

 - `c`     -- copy
"""
cp(p::AbstractVector{<:Real}) = copy(p)

"""
    cpv(pv::AbstractVector{<:AbstractVector{<:Real}})

Copy a position/velocity vector.

# Input

 - `pv`     -- pv-vector to be copied

# Output

 - `c`      -- copy
"""
cpv(pv::AbstractVector{<:AbstractVector{<:Real}}) = deepcopy(pv)

"""
    cr(r::AbstractMatrix{<:Real})

Copy an r-matrix.

# Input

 - `r`     -- r-matrix to be copied

# Output

 - `c`     -- copy
"""
cr(r::AbstractMatrix{<:Real}) = copy(r)

"""
    p2pv(p::AbstractVector{<:Real})

Extend a p-vector to a pv-vector by appending a zero velocity.

# Input

 - `p`     -- p-vector

# Output

 - `pv`    -- pv-vector
"""
function p2pv(p::AbstractVector{<:Real})
    zerot = zero(eltype(p))
    SVector(p, SVector(zerot, zerot, zerot))
end

"""
    pv2p(pv::AbstractVector{<:AbstractVector{<:Real}})

Discard velocity component of a pv-vector.

# Input

 - `pv`    -- pv-vector

# Output

 - `p`     -- p-vector
"""
pv2p(pv::AbstractVector{<:AbstractVector{<:Real}}) = pv[1]

#### Vector - Matrix / Initialization

"""
    ir()

Initialize an r-matrix to the identity matrix.

# Output

 - `r`     -- r-matrix
"""
ir() = SMatrix{3, 3}(1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, 1.0)

"""
    zp()

Zero a p-vector.

# Output
 - `p`     -- zero p-vector
"""
zp() = SVector{3}(zeros(3))

"""
    zpv()

Zero a pv-vector.

# Output

 - `pv`    -- zero pv-vector
"""
zpv() = SVector(zp(), zp())

"""
    zr()

Initialize an r-matrix to the null matrix.

# Output
 - `r`     -- r-matrix
"""
zr() = SMatrix{3, 3}(0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0)

#### Vector - Matrix / Matrix Operations

"""
    rxr(a::AbstractMatrix{<:Real}, b::AbstractMatrix{<:Real})

Multiply two r-matrices.

# Input

 - `a`     -- first r-matrix
 - `b`     -- second r-matrix

# Output

 - `atb`   -- a * b

# Note

1) It is permissible to re-use the same array for any of the
   arguments.
"""
rxr(a::V, b::V) where V<:AbstractMatrix{<:Real} = a*b

"""
    tr(r::AbstractMatrix{<:Real})

Transpose an r-matrix.

# Input

 - `r`     -- r-matrix

# Output

 - `rt`    -- transpose

# Note

1) It is permissible for r and rt to be the same array.
"""
tr(r::AbstractMatrix{<:Real}) = r'

#### Vector - Matrix / Matrix-Vector Products

"""
    rxp(r::AbstractMatrix{<:Real}, p::AbstractVector{<:Real})

Multiply a p-vector by an r-matrix.

# Input

 - `r`     -- r-matrix
 - `p`     -- p-vector

# Output

 - `rp`    -- r * p

# Note

1) It is permissible for p and rp to be the same array.
"""
rxp(r::AbstractMatrix{<:Real}, p::AbstractVector{<:Real}) = r*p

"""
    rxpv(r::AbstractMatrix{<:Real}, pv::AbstractVector{<:AbstractVector{<:Real}})

Multiply a pv-vector by an r-matrix.

# Input

 - `r`     -- r-matrix
 - `pv`    -- pv-vector

# Output

 - `rpv`   -- r * pv

# Note

1) The algorithm is for the simple case where the r-matrix r is not a
   function of time.  The case where r is a function of time leads to
   an additional velocity component equal to the product of the
   derivative of r and the position vector.

2) It is permissible for pv and rpv to be the same array.
"""
rxpv(r::AbstractMatrix{<:Real},
    pv::AbstractVector{<:AbstractVector{<:Real}}) = SVector{2}(r*pv[1], r*pv[2])

"""
    trxp(r::AbstractMatrix{<:Real}, p::AbstractVector{<:Real})

Multiply a p-vector by the transpose of an r-matrix.

# Input

 - `r`     -- r-matrix
 - `p`     -- p-vector

# Output

 - `trp`   -- r^T * p

# Note

1) It is permissible for p and trp to be the same array.
"""
trxp(r::AbstractMatrix{<:Real}, p::AbstractVector{<:Real}) = r'*p

"""
    trxpv(r::AbstractMatrix{<:Real}, pv::AbstractVector{<:AbstractVector{<:Real}})

Multiply a pv-vector by the transpose of an r-matrix.

# Input

 - `r`     -- r-matrix
 - `pv`    -- pv-vector

# Output

 - `trpv`  -- r^T * pv

# Note

1) The algorithm is for the simple case where the r-matrix r is not a
   function of time.  The case where r is a function of time leads to
   an additional velocity component equal to the product of the
   derivative of the transpose of r and the position vector.

2) It is permissible for pv and rpv to be the same array.
"""
trxpv(r::AbstractMatrix{<:Real},
    pv::AbstractVector{<:AbstractVector{<:Real}}) = SVector{2}(r'*pv[1], r'*pv[2])

#### Vector - Matrix / Rotation Vectors

"""
    rm2v(r::AbstractMatrix{<:Real})

Express an r-matrix as an r-vector.

# Input

 - `r`     -- rotation matrix

# Output

 - `w`     -- rotation vector (Note 1)

# Note

1) A rotation matrix describes a rotation through some angle about
   some arbitrary axis called the Euler axis.  The "rotation vector"
   returned by this function has the same direction as the Euler axis,
   and its magnitude is the angle in radians.  (The magnitude and
   direction can be separated by means of the function eraPn.)

2) If r is null, so is the result.  If r is not a rotation matrix the
   result is undefined; r must be proper (i.e. have a positive
   determinant) and real orthogonal (inverse = transpose).

3) The reference frame rotates clockwise as seen looking along the
   rotation vector from the origin.
"""
function rm2v(r::AbstractMatrix{<:Real})
    x, y, z = r[2,3] - r[3,2], r[3,1] - r[1,3], r[1,2] - r[2,1]
    s2, c2 = norm((x, y, z)), r[1,1] + r[2,2] + r[3,3] - 1
    zerot = zero(eltype(r))
    s2 > 0 ? SVector{3}(x, y, z)*atan(s2, c2)/s2 : SVector{3}(zerot, zerot, zerot)
end

"""
    rv2m(w::AbstractVector{<:Real})

Form the r-matrix corresponding to a given r-vector.

# Input

 - `w`     -- rotation vector (Note 1)

# Output

 - `r`     -- rotation matrix

# Note

1) A rotation matrix describes a rotation through some angle about
   some arbitrary axis called the Euler axis.  The "rotation vector"
   supplied to This function has the same direction as the Euler axis,
   and its magnitude is the angle in radians.

2) If w is null, the identity matrix is returned.

3) The reference frame rotates clockwise as seen looking along the
   rotation vector from the origin.
"""
function rv2m(w::AbstractVector{<:Real})
    #  Euler angle (magnitude of rotation vector)
    ϕ = norm(w)
    #  Euler axis (direction of rotation vector), perhaps null
    k = -(ϕ > 0 ? w/ϕ : w)
    I + sin(ϕ)*vec2mat(k) + (1-cos(ϕ))*vec2mat(k)*vec2mat(k)
end

#### Vector - Matrix / Separation and Angle

"""
    pap(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})

Position-angle from two p-vectors.

# Input

 - `a`     -- direction of reference point
 - `b`     -- direction of point whose PA is required

# Output

 - `θ`     -- position angle of b with respect to a (radians)

# Note

1) The result is the position angle, in radians, of direction b with
   respect to direction a.  It is in the range -pi to +pi.  The sense
   is such that if b is a small distance "north" of a the position
   angle is approximately zero, and if b is a small distance "east" of
   a the position angle is approximately +pi/2.

2) The vectors a and b need not be of unit length.

3) Zero is returned if the two directions are the same or if either
   vector is null.

4) If vector a is at a pole, the result is ill-defined.
"""
function pap(a::V, b::V) where V<:AbstractVector{<:Real}
    if norm(a) == 0 || norm(b) == 0
        θ = zero(eltype(a))
    else
        #  The north axis tangent from a (arbitrary length)
        η = SVector{3}(-a[1]*a[3], -a[2]*a[3], sum(a[1:2].^2))
        #  The east axis tanget from a (same length)
        ξ = vec2mat(η)*a/norm(a)
        # Resolve into components along the north and east axes
        θ = (b.-a)'*ξ == 0 && (b.-a)'*η == 0 ? zero(eltype(a)) :
            atan((b.-a)'*ξ, (b.-a)'*η)
    end
    θ
end

"""
    pas(λa::Real, ϕa::Real, λb::Real, ϕb::Real)

Position-angle from spherical coordinates.

# Input

 - `λa`    -- longitude of point A (e.g. RA) in radians
 - `ϕa`    -- latitude of point A (e.g. Dec) in radians
 - `λb`    -- longitude of point B
 - `ϕb`    -- latitude of point B

# Output

 - `θ`     -- position angle of B with respect to A

# Note

1) The result is the bearing (position angle), in radians, of point B
   with respect to point A.  It is in the range -pi to +pi.  The sense
   is such that if B is a small distance "east" of point A, the
   bearing is approximately +pi/2.

2) Zero is returned if the two points are coincident.
"""
function pas(λa::F, ϕa::F, λb::F, ϕb::F) where F<:Real
    x = sin(ϕb)*cos(ϕa) - cos(ϕb)*sin(ϕa)*cos(λb - λa)
    y = sin(λb - λa)*cos(ϕb)
    x != 0 || y != 0 ? atan(y, x) : 0.0
end

"""
    sepp(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})

Angular separation between two p-vectors.

# Input

 - `a`     -- first p-vector (not necessarily unit length)
 - `b`     -- second p-vector (not necessarily unit length)

# Output

 - `θ`     -- angular separation (radians, always positive)

# Note

1) If either vector is null, a zero result is returned.

2) The angular separation is most simply formulated in terms of scalar
   product.  However, this gives poor accuracy for angles near zero
   and pi.  The present algorithm uses both cross product and dot
   product, to deliver full accuracy whatever the size of the angle.
"""
function sepp(a::V, b::V) where V<:AbstractVector{<:Real}
    #  Sine of angle between the vectors, multiplied by the two moduli
    #  Cosine of the angle, multiplied by the two moduli
    cosθ, sinθ = sum(a.*b), norm(vec2mat(a)*b)
    sinθ != 0 || cosθ != 0 ? atan(sinθ, cosθ) : 0.0
end

"""
    seps(λa::Real, ϕa::Real, λb::Real, ϕb::Real)

Angular separation between two sets of spherical coordinates.

# Input

 - `λa`     -- first longitude (radians)
 - `ϕa`     -- first latitude (radians)
 - `λb`     -- second longitude (radians)
 - `ϕb`     -- second latitude (radians)

# Output

 - `θ`      -- angular separation (radians)

"""
function seps(λa::F, ϕa::F, λb::F, ϕb::F) where F<:Real
    #=
    #  Spherical to Cartesian
    a = [cos(λa)*cos(ϕa), sin(λa)*cos(ϕa), sin(ϕa)]
    b = [cos(λb)*cos(ϕb), sin(λb)*cos(ϕb), sin(ϕb)]
    #  Sine of angle between the vectors, multiplied by the two moduli
    sinθ = norm([0.0 -a[3] a[2]; a[3] 0.0 -a[1]; -a[2] a[1] 0.0]*b)
    #  Cosine of the angle, multiplied by the two moduli
    cosθ = sum(a.*b)
    sinθ != 0 || cosθ != 0 ? atan(sinθ, cosθ) : 0.0
    =#
    @inline sepp(s2c(λa, ϕa), s2c(λb, ϕb))
end

#### Vector - Matrix / Spherical-Cartesian

"""
    c2s(pos::AbstractVector{<:Real})

P-vector to spherical coordinates.

# Input

 - `p`     -- p-vector

# Output

 - `θ`     -- longitude angle (radians)
 - `ϕ`     -- latitude angle (radians)

# Note

1) The vector p can have any magnitude; only its direction is used.

2) If p is null, zero θ and ϕ are returned.

3) At either pole, zero θ is returned.
"""
function c2s(pos::AbstractVector{<:Real})
    zerot = zero(eltype(pos))
    NamedTuple{(:θ, :ϕ)}
    ((pos[1]^2 + pos[2]^2) == zerot ? zerot : atan(pos[2], pos[1]),
      pos[3] == zerot ? zerot : atan(pos[3], sqrt(pos[1]^2 + pos[2]^2)))
end

"""
    p2s(pos::AbstractVector{<:Real})

P-vector to spherical polar coordinates.

# Input

 - `p`     -- p-vector

# Output

 - `θ`     -- longitude angle (radians)
 - `ϕ`     -- latitude angle (radians)
 - `r`     -- radial distance

# Note

1) If P is null, zero θ, ϕ and r are returned.

2) At either pole, zero θ is returned.
"""
function p2s(pos::AbstractVector{<:Real})
    @inline NamedTuple{(:θ, :ϕ, :r)}((c2s(pos)..., norm(pos)))
end

"""
    pv2s(pv::AbstractVector{<:AbstractVector{<:Real}})

Convert position/velocity from Cartesian to spherical coordinates.

# Input

 - `posvel` -- position-velocity-vector

# Output

 - `θ`     -- longitude angle (radians)
 - `ϕ`     -- latitude angle (radians)
 - `r`     -- radial distance
 - `dθ`    -- rate of change of θ
 - `dϕ`    -- rate of change of ϕ
 - `dr`    -- rate of change of r

# Note

1) If the position part of pv is null, theta, ϕ, td and pd are
   indeterminate.  This is handled by extrapolating the position
   through unit time by using the velocity part of pv.  This moves the
   origin without changing the direction of the velocity component.
   If the position and velocity components of pv are both null, zeroes
   are returned for all six results.

2) If the position is a pole, theta, td and pd are indeterminate.  In
   such cases zeroes are returned for all three.
"""
function pv2s(pv::V) where V<:AbstractVector{<:AbstractVector{<:Real}}
    
    zerot = zero(eltype(pv[1]))
    x,  y,  z  = pv[norm(pv[1]) == zerot ? 2 : 1]
    dx, dy, dz = pv[2]

    if norm((x, y)) != zerot
        θ  = atan(y, x)
        ϕ  = atan(z, norm((x, y)))
        dθ = (x*dy - y*dx) / (x*x+y*y)
        dϕ = (dz*(x*x+y*y) - z*(x*dx+y*dy)) / (sum((x, y, z).^2)*norm((x, y)))
    else
        θ,  ϕ  = zerot, (z != zerot) ? atan(z, norm((x, y))) : zerot
        dθ, dϕ = zerot, zerot
    end
    r  = norm(pv[1])
    dr = norm((x, y, z)) != zerot ? (x*dx+y*dy+z*dz)/norm((x, y, z)) : zerot
    (; θ = θ, ϕ = ϕ, r = r, δθ = dθ, δϕ = dϕ, δr = dr)
end

"""
    s2c(θ::Real, ϕ::Real)

Convert spherical coordinates to Cartesian.

# Input

 - `θ`     -- longitude angle (radians)
 - `ϕ`     -- latitude angle (radians)

# Output

 - `c`     -- direction cosines

"""
function s2c(θ::F, ϕ::F) where F<:Real
    MVector(cos(θ)*cos(ϕ), sin(θ)*cos(ϕ), sin(ϕ))
end

"""
    s2p(θ::Real, ϕ::Real, r::Real)

Convert spherical polar coordinates to p-vector.

# Input

 - `θ`     -- longitude angle (radians)
 - `ϕ`     -- latitude angle (radians)
 - `r`     -- radial distance

# Output

 - `p`     -- Cartesian coordinates
"""
s2p(θ::F, ϕ::F, r::F) where F<:Real = r*s2c(θ, ϕ)

"""
    s2pv(θ::Real, ϕ::Real, r::Real, dθ::Real, dϕ::Real,
         dr::Real)

Convert position/velocity from spherical to Cartesian coordinates.

# Input

 - `θ`     -- longitude angle (radians)
 - `ϕ`     -- latitude angle (radians)
 - `r`     -- radial distance
 - `dθ`    -- rate of change of θ
 - `dϕ`    -- rate of change of ϕ
 - `dr`    -- rate of change of r

# Output

 - `pv`    -- pv-vector
"""
function s2pv(θ::F, ϕ::F, r::F, dθ::F, dϕ::F, dr::F) where F<:Real
    MVector(MVector(r*cos(θ)*cos(ϕ), r*sin(θ)*cos(ϕ), r*sin(ϕ)),
     MVector(-r*dθ*sin(θ)*cos(ϕ) - cos(θ)*(r*dϕ*sin(ϕ) - dr*cos(ϕ)),
      r*dθ*cos(θ)*cos(ϕ) - sin(θ)*(r*dϕ*sin(ϕ) - dr*cos(ϕ)),
      r*dϕ*cos(ϕ) + dr*sin(ϕ)))
end

#### Vector - Matrix / Vector Operations

"""
    pdp(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})

P-vector inner (=scalar=dot) product.

# Input

 - `a`     -- first p-vector
 - `b`     -- second p-vector

# Output

 - `r`     -- a . b
"""
pdp(a::V, b::V) where V<:AbstractVector{<:Real} = sum(a.*b)

"""
    pm(p::AbstractVector{<:Real}) = norm(p)

Modulus of p-vector.

# Input

 - `p`     -- p-vector

# Output

 - `r`     -- modulus
"""
pm(p::AbstractVector{<:Real}) = norm(p)

"""
    pmp(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})

P-vector subtraction.

# Input

 - `a`     -- first p-vector
 - `b`     -- second p-vector

# Output

 - `amb`   -- a - b
"""
pmp(a::V, b::V) where V<:AbstractVector{<:Real} = a .- b

"""
    pn(p::AbstractVector{<:Real})

Convert a p-vector into modulus and unit vector.

# Input

 - `p`     -- p-vector

# Output

 - `r`     -- modulus
 - `u`     -- unit vector

# Note

1) If p is null, the result is null.  Otherwise the result is a unit
   vector.
"""
function pn(p::AbstractVector{<:Real})
    zerot = zero(eltype(p))
    NamedTuple{(:modulus, :unit)}
    (norm(p) == 0 ? (zerot, SVector{3}(zerot, zerot, zerot)) :
        (norm(p), p./norm(p)))
end

"""
    ppp(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})

P-vector addition.

# Input

 - `a`     -- first p-vector
 - `b`     -- second p-vector

# Output

 - `apb`   -- a + b
"""
ppp(a::V, b::V) where V<:AbstractVector{<:Real} = a.+b

"""
    ppsp(a::AbstractVector{<:Real}, s::Real, b::AbstractVector{<:Real})

P-vector plus scaled p-vector.

# Input

 - `a`     -- first p-vector
 - `s`     -- scalar (multiplier for b)
 - `b`     -- second p-vector

# Output

 - `apsb`  -- a + s*b
"""
ppsp(a::V, s::Real, b::V) where V<:AbstractVector{<:Real} = a + s*b

"""
    pvdpv(a::AbstractVector{<:AbstractVector{<:Real}}, b::AbstractVector{<:AbstractVector{<:Real}})

Inner (=scalar=dot) product of two pv-vectors.

# Input

 - `a`     -- first pv-vector
 - `b`     -- second pv-vector

# Output

 - `adb`   -- a . b (see note)

# Note

1) If the position and velocity components of the two pv-vectors are (
   ap, av ) and ( bp, bv ), the result, a . b, is the pair of numbers
   ( ap . bp , ap . bv + av . bp ).  The two numbers are the
   dot-product of the two p-vectors and its derivative.
"""
function pvdpv(a::V, b::V) where
    V<:AbstractVector{<:AbstractVector{<:Real}}
   SVector{2}(sum(a[1].*b[1]), sum(a[1].*b[2] .+ a[2].*b[1]))
end

"""
    pvm(pv::AbstractVector{<:AbstractVector{<:Real}})

Modulus of pv-vector.

# Input

 - `pv`    -- pv-vector

# Output

 - `r`     -- modulus of position component
 - `s`     -- modulus of velocity component
"""
function pvm(pv::V) where
    V<:AbstractVector{<:AbstractVector{<:Real}}
    sqrt.(sum.(SVector{2}(pv[1].^2, pv[2].^2)))
end

"""
    pvmpv(a::AbstractVector{<:AbstractVector{<:Real}}, b::AbstractVector{<:AbstractVector{<:Real}})

Subtract one pv-vector from another.

# Input

 - `a`     -- first pv-vector
 - `b`     -- second pv-vector

# Output

 - `amb`   -- a - b
"""
function pvmpv(a::V, b::V) where
    V<:AbstractVector{<:AbstractVector{<:Real}}
    SVector{2}(a[1] .- b[1], a[2] .- b[2])
end

"""
    pvppv(a::AbstractVector{<:AbstractVector{<:Real}}, b::AbstractVector{<:AbstractVector{<:Real}})

Add one pv-vector to another.

# Input

 - `a`     -- first pv-vector
 - `b`     -- second pv-vector

# Output

 - `apb`   -- a + b
"""
function pvppv(a::V, b::V) where
    V<:AbstractVector{<:AbstractVector{<:Real}}
    SVector{2}(a[1] .+ b[1], a[2] .+ b[2])
end

"""
    pvu(dt::Real, pv::AbstractVector{<:AbstractVector{<:Real}})

Update a pv-vector.

# Input

 - `dt`    -- time interval
 - `pv`    -- pv-vector

# Output

 - `upv`   -- p updated, v unchanged

# Note

1) "Update" means "refer the position component of the vector to a new
   date dt time units from the existing date".

2) The time units of dt must match those of the velocity.
"""
function pvu(dt::Real, pv::V) where
    V<:AbstractVector{<:AbstractVector{<:Real}}
    SVector{2}(pv[1] .+ dt.*pv[2], pv[2])
end

"""
    pvup(dt::Real, pv::AbstractVector{<:AbstractVector{<:Real}})

Update a pv-vector, discarding the velocity component.

# Input

 - `dt`    -- time interval
 - `pv`    -- pv-vector

# Output

 - `p`     -- p-vector

# Note

1) "Update" means "refer the position component of the vector to a new
   date dt time units from the existing date".

2) The time units of dt must match those of the velocity.
"""
function pvup(dt::Real, pv::V) where
    V<:AbstractVector{<:AbstractVector{<:Real}}
    pv[1] .+ dt*pv[2]
end

"""
    pvxpv(a::AbstractVector{<:AbstractVector{<:Real}}, b::AbstractVector{<:AbstractVector{<:Real}})

Outer (=vector=cross) product of two pv-vectors.

# Input

 - `a`     -- first pv-vector
 - `b`     -- second pv-vector

# Output

 - `axb`   -- a x b

# Note

1) If the position and velocity components of the two pv-vectors are (
   ap, av ) and ( bp, bv ), the result, a x b, is the pair of vectors
   ( ap x bp, ap x bv + av x bp ).  The two vectors are the
   cross-product of the two p-vectors and its derivative.
"""
function pvxpv(a::V, b::V) where
    V<:AbstractVector{<:AbstractVector{<:Real}}
    SVector{2}(vec2mat(a[1])*b[1], vec2mat(a[1])*b[2] .+ vec2mat(a[2])*b[1])
end

"""
    pxp(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})

P-vector outer (=vector=cross) product.

# Input

 - `a`     -- first p-vector
 - `b`     -- second p-vector

# Output

 - `axb`   -- a x b
"""
function pxp(a::A, b::B) where {A<:AbstractVector{<:Real},
    B<:AbstractVector{<:Real}}
    vec2mat(a)*b
end

"""
    s2xpv(s1::Real, s2::Real, pv::AbstractVector{<:AbstractVector{<:Real}})

Multiply a pv-vector by two scalars.

# Input

 - `s1`    -- scalar to multiply position component by
 - `s2`    -- scalar to multiply velocity component by
 - `pv`    -- pv-vector

# Output

 - `spv`   -- pv-vector: p scaled by s1, v scaled by s2
"""
function s2xpv(s1::F, s2::F, pv::V) where
    {F<:Real, V<:AbstractVector{<:AbstractVector{<:Real}}}
    SVector{2}(s1*pv[1], s2*pv[2])
end

"""
    sxp(s::Real, p::AbstractVector{<:Real})

Multiply a p-vector by a scalar.

# Input

 - `s`     -- scalar
 - `p`     -- p-vector

# Output

 - `sp`    -- s * p
"""
sxp(s::Real, p::AbstractVector{<:Real}) = s*p

"""
    sxpv(s::Real, pv::AbstractVector{<:AbstractVector{<:Real}})

Multiply a pv-vector by a scalar.

# Input

 - `s`     -- scalar
 - `pv`    -- pv-vector

# Output

 - `spv`   -- s * pv
"""
function sxpv(s::F, pv::V) where
    {F<:Real, V<:AbstractVector{<:AbstractVector{<:Real}}}
    SVector{2}(s*pv[1], s*pv[2])
end
