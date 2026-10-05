# GyreInABox

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://aria-verify.github.io/GyreInABox.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://aria-verify.github.io/GyreInABox.jl/dev/)
[![Build Status](https://github.com/aria-verify/GyreInABox.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/aria-verify/GyreInABox.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/aria-verify/GyreInABox.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/aria-verify/GyreInABox.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

[Oceananigans](https://github.com/CliMA/Oceananigans.jl) based idealized models of ocean gyres in a bounded domain.

Currently two configurations are supported
- A model of wind and buoyancy forced ocean gyre adapted from the MITgcm [baroclinic ocean gyre example from documentation](https://mitgcm.readthedocs.io/en/latest/examples/baroclinic_gyre/baroclinic_gyre.html).
- An idealized model of the subpolar gyre with surface wind, temperature and salinity forcing and simplified topography with a northern basin separated from a southern open ocean region, based on the model described in [Spall (2011)](https://doi.org/10.1175/2011JCLI4130.1) and [Spall (2012)](https://doi.org/10.1175/JPO-D-11-0230.1).

A simple shared model interface is defined, with model-agnostic simulation driver, output handling and plotting built on top of this.
