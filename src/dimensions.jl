
abstract type AbstractSpatialDimension end
abstract type AbstractFixedSpatialDimension <: AbstractSpatialDimension end

struct FreeSpatialDimension <: AbstractSpatialDimension end
struct ReducedSpatialDimension <: AbstractFixedSpatialDimension end

struct SlicedSpatialDimension{T} <: AbstractFixedSpatialDimension
    coordinate::T
end

struct SpatialDimensions{X,Y,Z}
    x::X
    y::Y
    z::Z
end

function update(original::SpatialDimensions; x=nothing, y=nothing, z=nothing)
    SpatialDimensions(
        isnothing(x) ? original.x : x,
        isnothing(y) ? original.y : y,
        isnothing(z) ? original.z : z,
    )
end

index(::AbstractSpatialDimension, grid, nodes, N) = Colon()

function index(dimension::SlicedSpatialDimension, grid, nodes, N)
    clamp(searchsortedfirst(nodes(grid, Face()), dimension.coordinate), 1:N)
end

function indices(dimensions::SpatialDimensions, grid::RectilinearGrid)
    (
        index(dimensions.x, grid, xnodes, grid.Nx),
        index(dimensions.y, grid, ynodes, grid.Ny),
        index(dimensions.z, grid, znodes, grid.Nz),
    )
end

function indices(dimensions::SpatialDimensions, grid::LatitudeLongitudeGrid)
    (
        index(dimensions.x, grid, λnodes, grid.Nx),
        index(dimensions.y, grid, φnodes, grid.Ny),
        index(dimensions.z, grid, znodes, grid.Nz),
    )
end

function indices(dimensions::SpatialDimensions, grid::ImmersedBoundaryGrid)
    indices(dimensions, grid.underlying_grid)
end

const FixedX = SpatialDimensions{
    <:AbstractFixedSpatialDimension,<:FreeSpatialDimension,<:FreeSpatialDimension
}
const FixedY = SpatialDimensions{
    <:FreeSpatialDimension,<:AbstractFixedSpatialDimension,<:FreeSpatialDimension
}
const FixedZ = SpatialDimensions{
    <:FreeSpatialDimension,<:FreeSpatialDimension,<:AbstractFixedSpatialDimension
}
const FixedYZ = SpatialDimensions{
    <:FreeSpatialDimension,<:AbstractFixedSpatialDimension,<:AbstractFixedSpatialDimension
}
const FixedXZ = SpatialDimensions{
    <:AbstractFixedSpatialDimension,<:FreeSpatialDimension,<:AbstractFixedSpatialDimension
}
const FixedXY = SpatialDimensions{
    <:AbstractFixedSpatialDimension,<:AbstractFixedSpatialDimension,<:FreeSpatialDimension
}
const FixedXYZ = SpatialDimensions{
    <:AbstractFixedSpatialDimension,
    <:AbstractFixedSpatialDimension,
    <:AbstractFixedSpatialDimension,
}

const FixedYOrZ = Union{FixedY,FixedZ}
const FixedZOrXZ = Union{FixedZ,FixedXZ}
const FixedXOrYOrXY = Union{FixedX,FixedY,FixedXY}

const TwoSpatialDimensions = Union{FixedX,FixedY,FixedZ}
const ZeroOrOneSpatialDimensions = Union{FixedYZ,FixedXZ,FixedXY,FixedXYZ}

"""
    $(FUNCTIONNAME)(dimensions, grid)

Label for horizontal axis for plots of field with dimensions `dimensions` on grid `grid`.
    
$(TYPEDSIGNATURES)
"""
function axis_xlabel end

axis_xlabel(::FixedX, ::LatitudeLongitudeGrid) = "Latitude ϕ / ᵒ"
axis_xlabel(::FixedX, ::RectilinearGrid) = "Northward coordinate y / m"
axis_xlabel(::FixedYOrZ, ::LatitudeLongitudeGrid) = "Longitude λ / ᵒ"
axis_xlabel(::FixedYOrZ, ::RectilinearGrid) = "Eastward coordinate x / m"
axis_xlabel(::ZeroOrOneSpatialDimensions, ::AbstractUnderlyingGrid) = "Time / years"

function axis_xlabel(dimensions::SpatialDimensions, grid::ImmersedBoundaryGrid)
    axis_xlabel(dimensions, grid.underlying_grid)
end

"""
    $(FUNCTIONNAME)(dimensions, grid)

Label for vertical axis for plots of field with dimensions `dimensions` on grid `grid`.
    
$(TYPEDSIGNATURES)
"""
function axis_ylabel end

axis_ylabel(::FixedXOrYOrXY, ::AbstractUnderlyingGrid) = "Depth z / m"
axis_ylabel(::FixedZOrXZ, ::LatitudeLongitudeGrid) = "Latitude ϕ / ᵒ"
axis_ylabel(::FixedZOrXZ, ::RectilinearGrid) = "Northward coordinate y / m"
axis_ylabel(::FixedYZ, ::LatitudeLongitudeGrid) = "Longitude λ / ᵒ"
axis_ylabel(::FixedYZ, ::RectilinearGrid) = "Eastward coordinate x / m"
axis_ylabel(::FixedXYZ, ::AbstractUnderlyingGrid) = ""

function axis_ylabel(dimensions::SpatialDimensions, grid::ImmersedBoundaryGrid)
    axis_ylabel(dimensions, grid.underlying_grid)
end

"""
    $(FUNCTIONNAME)(dimensions)

Axis aspect ratio for field with dimensions `dimensions`.
    
$(TYPEDSIGNATURES)
"""

axis_aspect_ratio(::SpatialDimensions) = nothing
axis_aspect_ratio(::FixedZ) = DataAspect()

"""
    $(FUNCTIONNAME)(dimensions, grid, times)

Horizontal (x) axis limits for field with dimensions `dimension` on grid `grid` and times `times`.
    
$(TYPEDSIGNATURES)
"""
function axis_xlimits end

function axis_xlimits(
    ::ZeroOrOneSpatialDimensions, ::AbstractUnderlyingGrid, times::AbstractVector
)
    extrema(times / 365day)
end

function axis_xlimits(::FixedYOrZ, grid::LatitudeLongitudeGrid, ::AbstractVector)
    extrema(λnodes(grid, Face()))
end

function axis_xlimits(::FixedYOrZ, grid::RectilinearGrid, ::AbstractVector)
    extrema(xnodes(grid, Face()))
end

function axis_xlimits(::FixedX, grid::LatitudeLongitudeGrid, ::AbstractVector)
    extrema(φnodes(grid, Face()))
end

function axis_xlimits(::FixedX, grid::RectilinearGrid, ::AbstractVector)
    extrema(ynodes(grid, Face()))
end

function axis_xlimits(
    dimensions::SpatialDimensions, grid::ImmersedBoundaryGrid, times::AbstractVector
)
    axis_xlimits(dimensions, grid.underlying_grid, times)
end

"""
    $(FUNCTIONNAME)(dimensions, grid, times)

Vertical (y) axis limits for field with dimensions `dimension` on grid `grid` and times `times`.
    
$(TYPEDSIGNATURES)
"""
function axis_ylimits end

axis_ylimits(::FixedXYZ, ::AbstractUnderlyingGrid, ::AbstractVector) = nothing

function axis_ylimits(::FixedZOrXZ, grid::LatitudeLongitudeGrid, ::AbstractVector)
    extrema(φnodes(grid, Face()))
end

function axis_ylimits(::FixedZOrXZ, grid::RectilinearGrid, ::AbstractVector)
    extrema(ynodes(grid, Face()))
end

function axis_ylimits(::FixedYZ, grid::LatitudeLongitudeGrid, ::AbstractVector)
    extrema(λnodes(grid, Face()))
end

function axis_ylimits(::FixedYZ, grid::RectilinearGrid, ::AbstractVector)
    extrema(xnodes(grid, Face()))
end

function axis_ylimits(::FixedXOrYOrXY, grid::AbstractUnderlyingGrid, ::AbstractVector)
    extrema(znodes(grid, Face()))
end

function axis_ylimits(
    dimensions::SpatialDimensions, grid::ImmersedBoundaryGrid, times::AbstractVector
)
    axis_ylimits(dimensions, grid.underlying_grid, times)
end
