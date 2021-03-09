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
        "Additional Docstrings" => "every.md"
    ]
)


@info "edit docs"

stuff=readlines("nist_stuff/html_stuff.txt")

function fix_html(f)
    lines=readlines(f)

    f2 = open(f, "w")
    for line in lines
        if occursin("</head>", line )
            line2 = replace(line, "</head>" => stuff[1]*"</head>")
            write(f2, line2*"\n")
        else
            write(f2, line*"\n")
        end

    end
    close(f2)

end
    
for d in readdir("build")
    if isdir("build/$d")
        f = "build/$d/index.html"
        if isfile(f)
            fix_html(f)
        end
        if !isfile("build/$d/nist-combined.css")
            cp("nist_stuff/nist-combined.css", "build/$d/nist-combined.css")
        end
    end

end

fix_html("build/index.html")
if !isfile("build/nist-combined.css")
    cp("nist_stuff/nist-combined.css", "build/nist-combined.css")
end


@info "Deploy docs ..."

deploydocs(
    repo = "github.com/kfgarrity/TightlyBound.jl.git",
    devbranch = "main"
)
