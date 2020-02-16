using Documenter
using SphereSurfaceHistogram

makedocs(
    sitename = "SphereSurfaceHistogram",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == true),
    modules = [SphereSurfaceHistogram],
    pages = [
        "SphereSurfaceHistogram.jl" => "index.md",
        "Details.md",
        "Plotting.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/ffreyer/SphereSurfaceHistogram.jl.git"
)
