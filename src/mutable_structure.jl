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
    transformation::Bool = false
    benders::Bool = false
    uc::Bool = false
    one_cut::Bool = false
    split_sp0::Bool = false
    alpha::Int = 5
    max_iter::Int = 500
end