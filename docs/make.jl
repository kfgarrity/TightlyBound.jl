push!(LOAD_PATH,"../src/")
using Documenter, TightlyBound

#makedocs(sitename="TightlyBound.jl Documentation")

@info "Making documentation..."
makedocs(
  sitename="TightlyBound.jl Documentation",
  authors = "Kevin F. Garrity",
  format = Documenter.HTML(
      assets = ["assets/favicon.ico"],
  ),
#        prettyurls = get(ENV, "CI", nothing) == "true",
#        canonical = "https://oxfordcontrol.github.io/COSMO.jl/stable/",
#        assets = ["assets/favicon.ico"; "assets/github_buttons.js"; "assets/custom.css"],
#        analytics = "UA-134239283-1",
#  ),
  pages = [
        "Home" => "index.md",
        "User Guide" => Any[
        "AA" =>  "A.md",
        "BB" => "B.md",
        ],
        "Index" => "theindex.md",
    ]
)

@info "Deploy docs ..."

deploydocs(
    repo = "github.com/kfgarrity/TightlyBound.jl.git",
    devbranch = "main"
)
