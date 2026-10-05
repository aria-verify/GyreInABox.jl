"""
Sigmoidal depth profile which smoothly varies from `f_top` at `z >> z_0` to `f_bottom` at
`z << z_0` with length scale for smooth transition `ℓ` and offset for midpoint of
smooth transition `z_0`.

$(SIGNATURES)
"""
function sigmoidal_depth_profile(z, z_0, ℓ, f_bottom, f_top)
    f_bottom + (f_top - f_bottom) * (1 + tanh((z - z_0) / ℓ)) / 2
end

"""
Compute hyperbolically spaced grid face coordinate at index `k` for coordinate
discretized into `size` cells from lower value `lower` to upper value `upper`
with stretching factor `stretching_factor`.

$(SIGNATURES)
"""
function hyperbolically_spaced_faces(k, size, lower, upper, stretching_factor)
    lower +
    (upper - lower) * tanh(stretching_factor * (k - 1) / size) / tanh(stretching_factor)
end

"""
Smooth piecewise defined step function that is zero for negative `d`, one for `d > 1`
and smoothly interpolates between these values for `0 ≤ d ≤ 1`.
"""
function smooth_step(d)
    if d < 0.0
        0
    elseif d < 1
        3d^2 - 2d^3
    else
        1.0
    end
end

"""
Piecewise linear function, with constant value `y_0` for `t < t_0` and `t >= t_0 + 2Δt`,
linearly ramping from `y_0` to `y_0 + Δy` for `t_0 <= t < t_0 + Δt` and linearly ramping
from `y_0 + Δy` to `y_0` for `t_0 + Δt <= t < t_0 + 2Δt`.
"""
function triangular_ramp(t, t_0, Δt, y_0, Δy)
    if t_0 <= t < t_0 + Δt
        y_0 + Δy * (t - t_0) / Δt
    elseif t_0 + Δt <= t < t_0 + 2Δt
        y_0 + Δy - Δy * (t - Δt - t_0) / Δt
    else
        y_0
    end
end

"""
$(TYPEDEF)

Mask for horizontal circular region with center at `(x_center, y_center)`
and radius `radius`.
"""
struct HorizontalCircularRegionMask{T}
    x_center::T
    y_center::T
    radius::T
end

function (mask::HorizontalCircularRegionMask)(i, j, k, grid, field)
    x, y, z = node(i, j, k, field)
    (x - mask.x_center)^2 + (y - mask.y_center)^2 < mask.radius^2
end

function Base.summary(mask::HorizontalCircularRegionMask)
    "circular_region_at_x_$(mask.x_center)_y_$(mask.y_center)_radius_$(mask.radius)"
end

"""
$(SIGNATURES)

Convert string `s` in CamelCase to snake_case.
"""
camel_to_snake_case(s::String) = join(lowercase.(split(s, r"(?<=[a-z])(?=[A-Z])")), "_")
