using TensorNetworkAD2
using Documenter

DocMeta.setdocmeta!(TensorNetworkAD2, :DocTestSetup, :(using TensorNetworkAD2); recursive=true)

makedocs(;
    modules=[TensorNetworkAD2],
    authors="YidaiZhang",
    sitename="TensorNetworkAD2.jl",
    format=Documenter.HTML(;
        canonical="https://YidaiZhang.github.io/TensorNetworkAD2.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/YidaiZhang/TensorNetworkAD2.jl",
    devbranch="main",
)
