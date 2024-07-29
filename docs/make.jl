using Documenter
using Literate
using CVsim
using QuantumOpticsBase

builddir = "build"

### Convert examples to markdown using Literate.jl #joinpath(@__DIR__,
examplesdir = "./docs/src/examples"
markdowndir = "./docs/src/examples_md"

if !isdir(markdowndir)
    @info "Creating markdown output directory at \"", markdowndir, "\""
    mkdir(markdowndir)
end

exmpl = readdir(examplesdir,join=true)
if !isempty(exmpl)
    @info "Converting examples to markdown"
    for file in exmpl
        if isfile(file) && endswith(file, ".jl")
            Literate.markdown(file, markdowndir; flavor = Literate.DocumenterFlavor())
        end
    end
    exmpl = "examples_md/".*readdir(markdowndir)
else
    exmpl = ""
end

### Build the documentation
pages = [
        "index.md",
        "Installation" => "installation.md",
        "CVsim Parameters" => "parameters.md",
        "State preparation" => [
            "General structure" => "preparation/general.md",
            "GKP states" => "preparation/GKP.md",
            "Cat states" => "preparation/Cat.md"
        ],
        "State stabilization" => [
          "General structure" => "stabilization/general.md",
          "GKP states" => "stabilization/GKP.md"
        ],
        "Wigner quasiprobability" => [
            "Single system state" => "wigner/wigner_single.md",
            "Joint system state" => "wigner/wigner_joint.md",
            ],
        "Examples" => exmpl,
    ]


makedocs(
    format = Documenter.HTML(
        assets = ["assets/custom.css"],
        prettyurls = true,
    ),
    modules = [CVsim],
    build = builddir,
    sitename = "CVsim.jl",
    pages = pages
    )

deploydocs(
    repo = "github.com/irojkov-ph/CVsim.jl.git",
    )
