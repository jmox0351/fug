using Documenter
using fug

makedocs(
    sitename = "fug",
    format = Documenter.HTML(),
    modules = [fug]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
