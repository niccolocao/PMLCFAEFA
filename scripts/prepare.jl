using Pkg

root = normpath(joinpath(@__DIR__, ".."))
sourcedir = joinpath(root, "julia", "source")
project = joinpath(root, "julia", "PMLRBridge")
targetdir = joinpath(project, "src", "model")
mkpath(targetdir)

for name in ("MixCFAEFA.jl", "Sgr.jl", "PMLCFAEFA.jl")
    source = read(joinpath(sourcedir, name), String)
    if name == "PMLCFAEFA.jl"
        dependency_import = "using MixCFAEFA, Sgr,"
        occursin(dependency_import, source) || error("Cannot find module imports")
        source = replace(source, dependency_import => "using ..MixCFAEFA, ..Sgr,"; count=1)
        source = replace(source, "Optim.only_fg!" => "Optim.NLSolversBase.only_fg!")

        restart_loop = "while code_error == 0.0\n"
        occursin(restart_loop, source) || error("Cannot find restart loop")
        source = replace(source, restart_loop => "while code_error == 0.0 && iter_err <= maxNumRestart\n"; count=1)
        ending = "    return vcat(res,err_0_count)"
        occursin(ending, source) || error("Cannot find return statement")
        source = replace(source, ending => "    code_error == 0.0 && error(\"Maximum number of restarts reached\")\n" * ending; count=1)
    end
    write(joinpath(targetdir, name), source)
end

Pkg.activate(project)
Pkg.add([
    "CategoricalArrays", "Clustering", "Combinatorics", "Distributions",
    "ForwardDiff", "JLD2", "Optim", "SortingAlgorithms", "SpecialFunctions",
    "StatsBase", "LinearAlgebra", "Random", "Statistics",
])
Pkg.instantiate()
println("Julia sources prepared: ", project)
