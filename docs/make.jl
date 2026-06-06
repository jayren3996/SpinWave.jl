using Documenter
using SpinWave

DocMeta.setdocmeta!(SpinWave, :DocTestSetup, :(using SpinWave); recursive=true)

makedocs(
    sitename = "SpinWave.jl",
    authors = "JayRen and contributors",
    modules = [SpinWave],
    checkdocs = :none,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://jayren3996.github.io/SpinWave.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "Manual" => [
            "Conventions" => "manual/conventions.md",
            "Models" => "manual/models.md",
            "Spectra" => "manual/spectra.md",
        ],
        "Examples" => [
            "Ferromagnetic Chain" => "examples/ferromagnetic-chain.md",
            "Antiferromagnetic Chain" => "examples/antiferromagnetic-chain.md",
            "Square-Lattice Antiferromagnet" => "examples/square-lattice-antiferromagnet.md",
            "Anisotropic Ferromagnetic Chain" => "examples/anisotropic-ferromagnetic-chain.md",
            "Next-Nearest-Neighbor Chain" => "examples/next-nearest-neighbor-chain.md",
        ],
        "API Reference" => "reference/api.md",
    ],
)

if get(ENV, "CI", "false") == "true"
    deploydocs(
        repo = "github.com/jayren3996/SpinWave.jl.git",
        devbranch = "master",
    )
end
