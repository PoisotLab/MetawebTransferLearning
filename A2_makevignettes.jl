import Literate

corefiles = sort(filter(f -> startswith(f, r"\d\d"), readdir()))

vignetteconfig = Dict(
    "repo_root_url" => "https://github.com/PoisotLab/MetawebTransferLearning",
    "codefence" => Pair("````julia", "````"),
    "flavor" => Literate.FranklinFlavor(),
    "credit" => false
)

for corefile in corefiles
    Literate.markdown(corefile, "vignettes"; config=vignetteconfig)
end

README = readlines("_README.md")
for f in readdir("vignettes/")
    push!(README, "\n")
    append!(README, readlines(joinpath("vignettes", f)))
end
open("README.md", "w") do readme
    for line in README
        println(readme, line)
    end 
end