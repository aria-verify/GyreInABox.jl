"""
Sigmoidal depth profile which smoothly varies from `f_top` at `z = 0` to `f_bottom` at
`z = -depth` with length scale for smooth transition `ℓ` and offset for midpoint of
smooth transition `z_0`.

$(SIGNATURES)
"""
function sigmoidal_depth_profile(z, z_0, ℓ, f_bottom, f_top)
    f_bottom + (f_top - f_bottom) * (1  + tanh((z - z_0) / ℓ)) / 2
end

"""
Initial temperature at latitude `φ` and depth `z` with parameters `p`.

$(SIGNATURES)

## Details

Computes a latitude and depth dependent initial temperature with thermocline like
profile which smoothly varies from a latitude dependent surface temperature to a
constant deep ocean temperature with a sigmoidal profile.
"""
function initial_temperature(φ, z, p)
    T_surface = reference_surface_temperature(φ, p)
    sigmoidal_depth_profile(z, p.z_thermocline, p.ℓ_thermocline, p.T_abyssal, T_surface)
end

"""
Initial salinity at latitude `φ` and depth `z` with parameters `p`.

$(SIGNATURES)

## Details

Computes a latitude and depth dependent initial salinity with halocline like profile
which smoothly varies from a latitude dependent surface salinity to a constant
deep ocean salinity with a sigmoidal profile.
"""
function initial_salinity(φ, z, p)
    S_surface = reference_surface_salinity(φ, p)
    sigmoidal_depth_profile(z, p.z_halocline, p.ℓ_halocline, p.S_0, S_surface)
end
