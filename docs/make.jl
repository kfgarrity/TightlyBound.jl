push!(LOAD_PATH,"../src/")
using Documenter, TightlyBound

#makedocs(sitename="TightlyBound.jl Documentation")


#        prettyurls = get(ENV, "CI", nothing) == "true",
#        canonical = "https://oxfordcontrol.github.io/COSMO.jl/stable/",
#        assets = ["assets/favicon.ico"; "assets/github_buttons.js"; "assets/custom.css"],
#        analytics = "UA-134239283-1",
#  ),


@info "Making documentation..."
makedocs(
    sitename="TightlyBound.jl Documentation",
    authors = "Kevin F. Garrity",
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
    ),
    pages = [
        "Home" => "index.md",
        "User Guide" => Any[
            "Running Calculations" =>  "ug_run.md",
            "Fit Coefficients" => "ug_fit.md",
        ],
        "Core User Interface" => Any[
            "Structs" => "structs.md",
            "Functions" => "core.md",
        ],
        "Every Docstring" => "every.md"
    ]
)


@info "Deploy docs ..."

deploydocs(
    repo = "github.com/kfgarrity/TightlyBound.jl.git",
    devbranch = "main"
)
