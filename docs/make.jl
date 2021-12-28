using Documenter
using fug

push!(LOAD_PATH,"../src/")
makedocs(
    sitename = "Fenske-Underwood-Gilliand Documation",
    pages = ["Index" => "index.md"],
    format = Documenter.HTML(prettyurls = false),
    modules = [fug]
)


# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://github.com/jmox0351/fug.git",
    devbranch = "gh-pages"
)
