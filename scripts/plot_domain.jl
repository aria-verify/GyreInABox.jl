using GyreInABox
using CairoMakie

parameters = SpallDGParameters()
fig = GyreInABox.plot_domain_and_forcing(parameters)
save("domain_and_forcing.png", fig)
