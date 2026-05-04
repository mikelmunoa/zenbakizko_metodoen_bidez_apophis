using Documenter
using GravitationSimulation

makedocs(
    sitename = "GravitationSimulation.jl",
    authors  = "Mikel Muñoa Illarregi",
    modules  = [GravitationSimulation],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical  = "https://mikelmunoa.github.io/zenbakizko_metodoen_bidez_apophis/stable",
        edit_link  = "main",
    ),
    pages = [
        "Home"          => "index.md",
        "Physics"       => "physics.md",
        "API Reference" => "api.md",
        "Tutorials"     => [
            "Apophis 2029 flyby" => "tutorials/apophis_2029.md",
        ],
    ],
    checkdocs = :exports,
)

deploydocs(
    repo         = "github.com/mikelmunoa/zenbakizko_metodoen_bidez_apophis.git",
    push_preview = true,
)
