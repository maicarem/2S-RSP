using Parameters

@with_kw mutable struct MainPar
    @deftype String
    """
    Class of parameters
    """
    uc_strat::Int = 4
    uc_strat_4_limit::Int = 2000 # nb of user cuts for strategy 4
    uc_tolerance::Float64 = 0.01
    time_limit::Int = 3600
    @assert time_limit >= 0 # time_limit = 0 means infinity
    nthreads::Int = 4 # number of threads used in JuMP
end