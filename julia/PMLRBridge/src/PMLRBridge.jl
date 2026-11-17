module PMLRBridge

include("model/MixCFAEFA.jl")
include("model/Sgr.jl")
include("model/PMLCFAEFA.jl")

const LAST_RESULT = Ref(Float64[])
const LAST_ERROR = Ref("")

Base.@ccallable function pmlr_fit(
    xp::Ptr{Int64}, covp::Ptr{Cdouble},
    n::Cint, p::Cint, q::Cint, K::Cint, c::Cint, ncov::Cint,
    tol::Cdouble, maxit::Cint, maxrest::Cint, replix::Cint,
    fixthresholds::Cdouble, latentmeans::Cdouble,
    lp::Ptr{Cdouble}, has_l::Cint,
)::Cint
    try
        LAST_RESULT[] = Float64[]
        LAST_ERROR[] = ""
        X = copy(unsafe_wrap(Array, xp, (Int(n), Int(p)); own=false))
        covariates = copy(unsafe_wrap(Array, covp, (Int(n), Int(ncov)); own=false))
        loading_structure = has_l == 1 ? copy(unsafe_wrap(Array, lp, (Int(p), Int(q)); own=false)) : false
        result = PMLCFAEFA.PMLcfaefa(
            X, Int(p), Int(q), Int(K), Int(c), covariates,
            Float64(tol), Int(maxit), Int(maxrest), Int(replix),
            Float64(fixthresholds), loading_structure, Float64(latentmeans), false,
        )
        LAST_RESULT[] = Float64.(result)
        return Cint(0)
    catch err
        LAST_ERROR[] = sprint(showerror, err, catch_backtrace())
        LAST_RESULT[] = Float64[]
        return Cint(1)
    end
end

Base.@ccallable function pmlr_result_length()::Cint
    return Cint(length(LAST_RESULT[]))
end

Base.@ccallable function pmlr_copy_result(dest::Ptr{Cdouble}, capacity::Cint)::Cint
    values = LAST_RESULT[]
    if capacity < length(values) || dest == C_NULL
        return Cint(1)
    end
    GC.@preserve values begin
        unsafe_copyto!(dest, pointer(values), length(values))
    end
    return Cint(0)
end

Base.@ccallable function pmlr_last_error(dest::Ptr{UInt8}, capacity::Cint)::Cint
    if capacity <= 0 || dest == C_NULL
        return Cint(0)
    end
    bytes = codeunits(LAST_ERROR[])
    count = min(length(bytes), Int(capacity) - 1)
    for i in 1:count
        unsafe_store!(dest, bytes[i], i)
    end
    unsafe_store!(dest, UInt8(0), count + 1)
    return Cint(count)
end

end
