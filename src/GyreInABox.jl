"""
Oceananigans based idealized models of ocean gyres in a bounded domain.

A simple shared model interface is defined, with model-agnostic simulation driver,
output handling and plotting built on top of this.
"""
module GyreInABox

using DocStringExtensions
using Oceananigans
using Oceananigans.DistributedComputations
using Oceananigans.Grids
using Oceananigans.Units
using SeawaterPolynomials: TEOS10EquationOfState
using CairoMakie
using Printf
using JLD2

@template FUNCTIONS = """
                      $(DOCSTRING)

                      $(TYPEDSIGNATURES)
                      """

@template METHODS = """
                    $(SIGNATURES)

                    $(DOCSTRING)
                    """

@template TYPES = """
                  $(TYPEDEF)

                  $(DOCSTRING)

                  $(TYPEDFIELDS)
                  """

export DoubleGyreParameters, Spall2011Parameters
export SimulationConfiguration
export EastwardVelocity, NorthwardVelocity, UpwardVelocity
export Salinity, Temperature
export SpecificKineticEnergy, SpecificTurbulentKineticEnergy
export EastwardBarotropicVelocity, NorthwardBarotropicVelocity, FreeSurfaceDisplacement
export MOCStreamFunction, BarotropicStreamFunction
export NorthwardHeatTransport
export FieldModelOutput, PointModelOutput
export horizontal_slice_output, x_depth_slice_output, y_depth_slice_output
export depth_averaged_output, horizontally_averaged_output, spatially_averaged_output
export free_surface_output, stream_functions_output
export moc_strength_at_y_output, northward_heat_transport_at_y_output
export setup_model, initialize!, setup_simulation
export AnimationPlotOutput, TemporalAveragePlotOutput, TimeSeriesPlotOutput
export run_simulation, plot_outputs

const gravitational_acceleration = 9.80665

include("variables.jl")
include("dimensions.jl")
include("processors.jl")
include("outputs.jl")
include("plots.jl")
include("models.jl")
include("simulations.jl")
include("utils.jl")
include("double_gyre.jl")
include("spall_2011.jl")

end
