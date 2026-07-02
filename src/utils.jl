"""
Sigmoidal depth profile which smoothly varies from `f_top` at `z = 0` to `f_bottom` at
`z = -depth` with length scale for smooth transition `ℓ` and offset for midpoint of
smooth transition `z_0`.

$(SIGNATURES)
"""
function sigmoidal_depth_profile(z, z_0, ℓ, f_bottom, f_top)
    f_bottom + (f_top - f_bottom) * (1 + tanh((z - z_0) / ℓ)) / 2
end

"""
Compute hyperbolically spaced grid face coordinate at index `k` for coordinate
discretized into `size` cells from lower value `lower` to upper value `upper`
with strething factor `stretching_factor`.

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
    if d < 0.
        0
    elseif d < 1
        3d^2 - 2d^3
    else
        1.
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
