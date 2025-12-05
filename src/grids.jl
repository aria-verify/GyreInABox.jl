"""
Compute hyperbolically spaced depth grid face coordinate at index `k` given model
configuration `configuration`.

$(SIGNATURES)
"""
function hyperbolically_spaced_faces(k::Int, configuration::GyreInABoxConfiguration)
    configuration.depth_interval[2] + (
        configuration.depth_interval[1] - configuration.depth_interval[2]
    ) * (
        1 - tanh(
            configuration.depth_stretching_factor * (k - 1) / configuration.grid_size[3]
        ) / tanh(configuration.depth_stretching_factor)
    )
end

"""
Set up ocean gyre model grid with configuration `configuration` optionally overriding
architecture in `configuration` with `architecture`.

$(SIGNATURES)
"""
function setup_grid(
    configuration::GyreInABoxConfiguration{T}; architecture = nothing
) where {T}
    LatitudeLongitudeGrid(
        isnothing(architecture) ? configuration.architecture : architecture,
        size=configuration.grid_size,
        longitude=configuration.longitude_interval,
        latitude=configuration.latitude_interval,
        z=k -> hyperbolically_spaced_faces(k, configuration),
        halo=configuration.halo_size,
        topology=(Bounded, Bounded, Bounded)
    )
end
