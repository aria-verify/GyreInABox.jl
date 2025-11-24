using GyreInABox
using Documenter

DocMeta.setdocmeta!(GyreInABox, :DocTestSetup, :(using GyreInABox); recursive=true)

makedocs(;
    modules=[GyreInABox],
    authors="Matt Graham and contributors",
    sitename="GyreInABox.jl",
    format=Documenter.HTML(;
        canonical="https://aria-verify.github.io/GyreInABox.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/aria-verify/GyreInABox.jl",
    devbranch="main",
)
