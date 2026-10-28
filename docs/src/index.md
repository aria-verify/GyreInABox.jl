```@meta
CurrentModule = GyreInABox
```

# GyreInABox

Documentation for [GyreInABox](https://github.com/aria-verify/GyreInABox.jl) which implements [Oceananigans](https://github.com/CliMA/Oceananigans.jl) based idealized models of ocean gyres in a bounded domain.

Currently two configurations are supported
- [A model of wind and buoyancy forced ocean gyre](models/double_gyre.md) adapted from the MITgcm [baroclinic ocean gyre example from documentation](https://mitgcm.readthedocs.io/en/latest/examples/baroclinic_gyre/baroclinic_gyre.html).
- An idealized model of the subpolar gyre with surface wind, temperature and salinity forcing and simplified topography with a northern basin separated from a southern open ocean region, based on the model described in [Spall (2011)](https://doi.org/10.1175/2011JCLI4130.1) and [Spall (2012)](https://doi.org/10.1175/JPO-D-11-0230.1).

A simple shared model interface is defined, with model-agnostic simulation driver, output handling and plotting built on top of this.

## API reference

Reference documentation for the package's public and internal interfaces are linked to below

```@contents
Pages = ["lib/public.md", "lib/internals.md"]
Depth = 1
```
