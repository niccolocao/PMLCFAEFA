using Pkg

root = normpath(joinpath(@__DIR__, ".."))
project = joinpath(root, "julia", "PMLRBridge")
isfile(joinpath(project, "src", "model", "PMLCFAEFA.jl")) || error("Run scripts/prepare.jl first")
Pkg.activate(project)
Pkg.instantiate()
Pkg.activate(; temp=true)
Pkg.add("PackageCompiler")
using PackageCompiler
Pkg.activate(project)

output = joinpath(root, "Rpkg", "inst", "julia")
mkpath(dirname(output))
PackageCompiler.create_library(
    project, output;
    lib_name="pmlr",
    force=true,
    header_files=[joinpath(root, "julia", "pmlr.h")],
    cpu_target="generic",
)
println("Compiled bundle ready in: ", output)
