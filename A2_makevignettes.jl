import Literate

corefiles = sort(filter(f -> startswith(f, r"\d\d"), readdir()))

vignetteconfig = Dict(
    "repo_root_url" => "https://github.com/PoisotLab/MetawebTransferLearning",
    "codefence" => Pair("````julia", "````"),
    "flavor" => Literate.FranklinFlavor()

)

for corefile in corefiles
    Literate.markdown(corefile, "vignettes"; config=vignetteconfig)
end