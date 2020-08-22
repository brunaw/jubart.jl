using jubart
using Documenter

makedocs(;
    modules=[jubart],
    authors="Bruna Wundervald, Mateus Maia",
    repo="https://github.com/brunaw/jubart.jl/blob/{commit}{path}#L{line}",
    sitename="jubart.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://brunaw.github.io/jubart.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/brunaw/jubart.jl",
)
