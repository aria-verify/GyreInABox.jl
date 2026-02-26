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

