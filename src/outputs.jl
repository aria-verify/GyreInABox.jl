const DEFAULT_SCHEDULE = TimeInterval(1day)

abstract type AbstractModelOutput{S} end

struct PointModelOutput{V,S,T} <: AbstractModelOutput{S}
    variables::V
    schedule::S
    points::Matrix{T}
end

"""
$(TYPEDEF)

Specification of model outputs to write out during simulation.

$(TYPEDFIELDS)
"""
struct FieldModelOutput{V,S,P} <: AbstractModelOutput{S}
    "Model variables to record as part of output"
    variables::V
    "Schedule to record outputs at"
    schedule::S
    "Optional spatial processor(s) to apply to fields associated with variables"
    processor::P
end

"""
$(SIGNATURES)

Horizontal (latitude-longitude or x-y) slice output.
    
## Details

Records horizontal slices through model fields at specified depth.
"""
function horizontal_slice_output(;
    depth, schedule=DEFAULT_SCHEDULE, variables=VELOCITY_AND_TRACER_VARIABLES
)
    FieldModelOutput(
        variables, schedule, SpatialSliceProcessor(; z=SlicedSpatialDimension(depth))
    )
end

"""
$(SIGNATURES)

Vertical (x-depth or longitude-depth) slice output.
    
## Details

Records vertical slices through model fields at specified northward y coordinate
(if on a rectilinear grid) or latitude coordinate (if on a latitude-longitude grid).
"""
function x_depth_slice_output(;
    y_or_latitude, schedule=DEFAULT_SCHEDULE, variables=VELOCITY_AND_TRACER_VARIABLES
)
    FieldModelOutput(
        variables,
        schedule,
        SpatialSliceProcessor(; y=SlicedSpatialDimension(y_or_latitude)),
    )
end

"""
$(SIGNATURES)

Vertical (y-depth or latitude-depth) slice output.
    
## Details

Records vertical slices through model fields at specified eastward x coordinate
(if on a rectilinear grid) or longitude coordinate (if on a latitude-longitude grid).
"""
function y_depth_slice_output(;
    x_or_longitude, schedule=DEFAULT_SCHEDULE, variables=VELOCITY_AND_TRACER_VARIABLES
)
    FieldModelOutput(
        variables,
        schedule,
        SpatialSliceProcessor(; x=SlicedSpatialDimension(x_or_longitude)),
    )
end

"""
$(SIGNATURES)

Free surface fields output.
    
## Details

Records two-dimensional free surface (displacement and barotropic velocity) fields.
"""
function free_surface_output(; schedule=DEFAULT_SCHEDULE)
    FieldModelOutput(
        (FreeSurfaceDisplacement(), BAROTROPIC_VELOCITY_VARIABLES...), schedule, nothing
    )
end

"""
$(SIGNATURES)

Depth averaged output.

## Details

Records horizontal fields corresponding to depth averaged model fields
computed over whole domain or a region specified by a binary mask.
"""
function depth_averaged_output(;
    schedule=DEFAULT_SCHEDULE, variables=VELOCITY_AND_TRACER_VARIABLES, mask=nothing
)
    FieldModelOutput(variables, schedule, SpatialAverageProcessor((3,), mask))
end

"""
$(SIGNATURES)

Horizontally averaged fields output
    
## Details

Records one-dimensional fields corresponding to horizontal averages of fields
at different depths, with horizontal averages computed over whole domain or a
region specified by a binary mask.
"""
function horizontally_averaged_output(;
    schedule=DEFAULT_SCHEDULE, variables=TRACER_VARIABLES, mask=nothing
)
    FieldModelOutput(variables, schedule, SpatialAverageProcessor((1, 2), mask))
end

"""
$(SIGNATURES)

Spatially averaged output.
    
## Details

Records scalar fields corresponding to average over whole domain or a region
specified by a binary mask.
"""
function spatially_averaged_output(;
    schedule=DEFAULT_SCHEDULE, variables=TRACER_VARIABLES, mask=nothing
)
    FieldModelOutput(variables, schedule, SpatialAverageProcessor((1, 2, 3), mask))
end

"""
$(SIGNATURES)

Meridional overturning circulation (MOC) and barotropic stream functions output.

## Details

Records two-dimensional fields corresponding to MOC and barotropic stream functions.

See [`MOCStreamFunction`](@ref) and [`BarotropicStreamFunction`](@ref) for more details
"""
function stream_functions_output(; schedule=DEFAULT_SCHEDULE)
    FieldModelOutput((MOCStreamFunction(), BarotropicStreamFunction()), schedule, nothing)
end

"""
$(SIGNATURES)

Meridional overturning circulation (MOC) strength output.
    
## Details

Records scalar time series corresponding to maximum over depth of stream function of 
meridional overturning circulation (MOC) at a given y / latitude coordinate.

The scalar recorded here corresponds to ``\\max_z \\Psi^M(\\varphi_r, z, t)``
(on a latitude-longitude grid) or ``\\max_z \\Psi^M(y_r, z, t)`` (on a rectilinear
grid) where ``\\varphi_r`` and ``y_r`` are the latitude and y coordinates the
strength is measured at respectively, where ``\\Psi^M`` is computed as described in
[`MOCStreamFunction`](@ref).
"""
function moc_strength_at_y_output(; y_or_latitude, schedule=DEFAULT_SCHEDULE)
    ModelOutput(
        (MOCStreamFunction(),),
        schedule,
        (
            SpatialMaximumProcessor((3,)),
            SpatialSliceProcessor(; y=SlicedSpatialDimension(y_or_latitude)),
        ),
    )
end

"""
$(SIGNATURES)

Northward heat transport output.
    
## Details

Records scalar fields corresponding to northward (meridional) heat transport at a
given y / latitude coordinate.

The northward heat transport is computed here as described in [`NorthwardHeatTransport`](@ref).
"""
function northward_heat_transport_at_y_output(; y_or_latitude, schedule=DEFAULT_SCHEDULE)
    FieldModelOutput(
        (NorthwardHeatTransport(),),
        schedule,
        SpatialSliceProcessor(; y=SlicedSpatialDimension(y_or_latitude)),
    )
end

"""
$(SIGNATURES)

Symbol label for output type `output` to use in naming output file and registering output writer.
"""
function label(output::FieldModelOutput)
    base_label =
        "variables_" *
        join((short_name(v) for v in output.variables), "_") *
        "_at_" *
        label(output.schedule)
    isnothing(output.processor) ? base_label : label(output.processor) * "_of_" * base_label
end

function label(schedule::Oceananigans.Utils.AbstractSchedule)
    replace(
        camel_to_snake_case(summary(schedule)),
        " " => "",
        "(" => "_",
        ")" => "",
        "," => "_",
        "=" => "_",
    )
end

"""
$(SIGNATURES)

Spatial grid indices output type `output` records fields at for grid `grid`.
"""
function indices(output::FieldModelOutput, grid)
    per_variable_indices = [
        indices(spatial_dimensions(variable, output.processor), grid) for
        variable in output.variables
    ]
    !allequal(per_variable_indices) && error(
        "Non-compatible indices computed for variables associated with output $(output)"
    )
    first(per_variable_indices)
end

indices(output::AbstractModelOutput, grid) = (:, :, :)

"""
$(SIGNATURES)

Named tuple of output fields deriving from those in `model` to 
record for output `model_output`.
"""
function fields(model_output::FieldModelOutput, model)
    NamedTuple(
        Symbol(short_name(variable)) =>
            process(model_output.processor, field(variable, model)) for
        variable in model_output.variables
    )
end

"""
$(SIGNATURES)

Time schedule to record output type at.
"""
schedule(output::AbstractModelOutput) = output.schedule

"""
$(SIGNATURES)

Filename to record outputs to.

## Details

For an output `output` a label computed using `label` function is appended on to
`stem`` and file extension is specified by `extension` added.
"""
function output_filename(stem::String, label::String, extension::String)
    "$(stem)_$(label).$(extension)"
end

function output_filename(stem::String, label::String, extension::Nothing)
    "$(stem)_$(label)"
end

function output_filename(
    stem::String, output::AbstractModelOutput, extension::Union{String,Nothing}=nothing
)
    output_filename(stem, label(output), extension)
end

"""
    $(FUNCTIONNAME)(output_writer_type)

File extension to use with `output_writer_type`.
    
$(TYPEDSIGNATURES)
"""
function extension end

extension(::Type{JLD2Writer}) = "jld2"
extension(::Type{NetCDFWriter}) = "nc"
extension(::Type{ZarrWriter}) = "zarr"
