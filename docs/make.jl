using GyreInABox
using Documenter

DocMeta.setdocmeta!(GyreInABox, :DocTestSetup, :(using GyreInABox); recursive = true)

makedocs(;
    modules = [GyreInABox],
    authors = "Matt Graham and contributors",
    sitename = "GyreInABox.jl",
    format = Documenter.HTML(;
        canonical = "https://aria-verify.github.io/GyreInABox.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = Any[
        "Home"=>"index.md",
        "Models"=>[
            "Wind and buoyancy forced double gyre"=>"models/double_gyre.md",
            "Wind, temperature and salinity forced gyre with basin bathymetry"=>"models/spall_2011.md",
        ],
        "API reference"=>["Public"=>"lib/public.md", "Internals"=>"lib/internals.md"],
    ],
    doctest = false,
)

deploydocs(; repo = "github.com/aria-verify/GyreInABox.jl", devbranch = "main")
