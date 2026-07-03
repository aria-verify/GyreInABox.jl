abstract type AbstractSpatialProcessor end

struct SpatialSliceProcessor{D} <: AbstractSpatialProcessor
    sliced_dimensions::D
end

function SpatialSliceProcessor(; kwargs...)
    SpatialSliceProcessor(NamedTuple(k => v for (k, v) in pairs(kwargs)))
end

function label(processor::SpatialSliceProcessor)
    "spatial_slice_at_" * join(
        (
            "$(name)_$(value.coordinate)" for
            (name, value) in pairs(processor.sliced_dimensions)
        ),
        "_",
    )
end

process(::Nothing, field_or_dimension) = field_or_dimension

process(::SpatialSliceProcessor, field::Field) = field

function process(processor::SpatialSliceProcessor, dimensions::SpatialDimensions)
    update(dimensions; (k => v for (k, v) in pairs(processor.sliced_dimensions))...)
end

abstract type AbstractSpatialReductionProcessor{D} <: AbstractSpatialProcessor end

struct SpatialAverageProcessor{D,M} <: AbstractSpatialReductionProcessor{D}
    dimensions::D
    mask::M
end

function SpatialAverageProcessor(dimensions::Tuple{Vararg{Int}})
    SpatialAverageProcessor(dimensions, nothing)
end

function label(processor::SpatialAverageProcessor)
    base_label = "spatial_average_over_" * join("xyz"[d] for d in processor.dimensions)
    isnothing(processor.mask) ? base_label : "$(base_label)_in_$(summary(processor.mask))"
end

function process(processor::SpatialAverageProcessor, field::Field)
    Field(Average(field; dims=processor.dimensions, condition=processor.mask))
end

struct SpatialMaximumProcessor{D} <: AbstractSpatialReductionProcessor{D}
    dimensions::D
end

function label(processor::SpatialMaximumProcessor)
    "spatial_maximum_over_" * join("xyz"[d] for d in processor.dimensions)
end

function process(processor::SpatialMaximumProcessor, field::Field)
    Field(Reduction(maximum!, field; dims=processor.dimensions))
end

function process(
    processor::AbstractSpatialReductionProcessor, dimensions::SpatialDimensions
)
    update(
        dimensions;
        (Symbol("xyz"[d]) => ReducedSpatialDimension() for d in processor.dimensions)...,
    )
end

function label(processor::Vector{<:AbstractSpatialProcessor})
    join((label(p) for p in processor), "_")
end

function process(processor::Vector{<:AbstractSpatialProcessor}, field::Field)
    for p in processor
        field = process(p, field)
    end
    field
end

function process(
    processor::Vector{<:AbstractSpatialProcessor}, dimensions::SpatialDimensions
)
    for p in processor
        dimensions = process(p, dimensions)
    end
    dimensions
end

function spatial_dimensions(variable::AbstractModelVariable, processor)
    process(processor, spatial_dimensions(variable))
end
