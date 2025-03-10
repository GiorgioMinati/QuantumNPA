function sym_add!(matrix, i, j, val)
    matrix[i, j] += val

    if i != j
        matrix[j, i] += val
    end

    return matrix
end



"""
Construct the NPA moment matrix.

The argument operators can in general be an array of arrays of operators
(blocks), e.g. [[Id], [A1 + B1], [A1 - B1]]. It can also be a simple array of
operators, in which case it is treated the same as an array containing a
single array of operators, e.g., [[Id, A1, A2]]). In either case the return
value is a dictionary with:

  * as keys: monomials obtained by multiplying operators in the same blocks
    together.

  * as values: block-diagonal sparse matrices with coefficients obtained
    from multiplying the input operators together.
"""
function npa_moment(operators::Vector{<:Union{Monomial,Polynomial}})
    N = length(operators)
    iops = collect(enumerate(operators))
    moment = Polynomial((N, N))

    for (i, x) in iops
        for (j, y) in iops[i:end]
            p = Polynomial(conj_min(conj(x)*y))

            for (c, m) in p
                if !hasmonomial(moment, m)
                    moment[m] = sym_add!(spzeros(N, N), i, j, c)
                else
                    sym_add!(moment[m], i, j, c)
                end
            end
        end
    end

    return moment
end

#function npa_moment(operators::Vector{Vector{T}} where T)
#    moment = npa_moment.(operators)
#    return blockdiag(moment, (sz) -> spzeros(Float64, sz))
#end

npa_moment(source, level) = npa_moment(npa_level(source, level))


"""
Generate the NPA relaxation for a given quantum optimisation problem (an
operator expr whose expectation we want to maximise with the expectation
values of the operator constraints set to zero).

The result is 
"""
function npa2sdp(expr, level; eq=[], ge=[])
    moment = npa_moment([expr, eq, ge], level)
    return npa2sdp(expr, moment, eq=eq, ge=ge)
end

function _reduce_problem(expr, moment::Polynomial; eq=[], ge=[])
    # Reduce constraints to canonical form
    expr = conj_min(expr)
    eq = linspace(map(conj_min, eq))
    ge = map(conj_min, ge)

    if haskey(eq, Id)
        @error "Contradiction Id = 0 in equality constraints."
    end

    # Reduce the objective expression, using constraints to eliminate
    # monomials.
    expr = reduce_expr(expr, eq)

    # Reduce moments using equality constraints.
    moment = reduce_expr(moment, eq)

    # Reduce inequality constraints then include them as inequalities along
    # with the original moment matrix.
    ge = reduce_exprs(ge, eq)
    
    return expr, moment, ge, eq
end

function npa2sdp(expr, moment::Polynomial; eq=[], ge=[])
    expr, moment, ge, _ = _reduce_problem(expr, moment, eq=eq, ge=ge) 
    return (expr, vcat([moment], ge))
end

if !@isdefined(default_solver)
    default_solver = SCS.Optimizer
end

function set_solver!(solver)
    global default_solver = solver
end

function set_verbosity!(model, verbose)
    if !isnothing(verbose)
        (!verbose ? set_silent : unset_silent)(model)
    end
end

function sparse_trprod(A::SparseMatrixCSC{Float64, Int64}, B::Symmetric{VariableRef, Matrix{VariableRef}}
)
    n = size(A)[1]
    if size(B)[2] != n
        @error "Matrix sizes do not match: A:$n, B:$(size(B)[2])"
    end
    Is, Js, Vs = findnz(A)

    exp = AffExpr(0)
    for k=1:length(Is)
        add_to_expression!(exp, Vs[k] * B[Js[k],Is[k]])
    end
    return exp
end

function sdp2jump(expr, ineqs;
                  goal=:maximise,
                  solver=nothing,
                  verbose=nothing)
    if goal in (:maximise, :maximize, :max)
        maximise = true
        s = 1
    elseif goal in (:minimise, :minimize, :min)
        maximise = false
        s = -1
    end
    
    model = !isnothing(solver) ? Model(solver) : Model()
    
    Zs = [@variable(model, [1:m, 1:n], PSD, set_string_name = false)
          for (i, (m, n)) in enumerate(size_as_pair.(ineqs))]
    model[:Zs] = Zs
    
    Ids = (ineq[Id] for ineq in ineqs)
    objective = (sum(LinearAlgebra.tr(s*m*z)
                     for (m, z) in zip(Ids, Zs))
                 + expr[Id])

    if maximise
        @objective(model, Min, objective)
    else
        @objective(model, Max, objective)
    end

    mons = collect(m for m in monomials(expr, ineqs) if !isidentity(m))

    for m in mons
        c = expr[m]
        Fs = (ineq[m] for ineq in ineqs)
        tr_term = @expression(model, sum(sparse_trprod(F, Z)
                                         for (F, Z) in zip(Fs, Zs)))
        
        @constraint(model, tr_term + s*c == 0, set_string_name = false)
    end

    #@expression(model, moncons[i=1:length(mons)], 
    #    sum(LinearAlgebra.tr(F*Z)
    #        for (F, Z) in zip((ineq[mons[i]] for ineq in ineqs), Zs)) 
    #    + expr[mons[i]]*s)
    #@constraint(model, moncons .== 0, set_string_name = false)

    set_verbosity!(model, verbose)

    return model
end

function npa2jump(expr, level_or_moments;
                  eq=[],
                  ge=[],
                  goal=:maximise,
                  solver=nothing,
                  verbose=nothing)
    (expr, moments) = npa2sdp(expr, level_or_moments, eq=eq, ge=ge)

    model = sdp2jump(expr, moments,
                     goal=goal,
                     solver=solver,
                     verbose=verbose)
    return model
end



function npa_opt(expr, level_or_moments;
                 eq=[],
                 ge=[],
                 goal=:maximise,
                 solver=default_solver,
                 verbose=false)
    model = npa2jump(expr, level_or_moments,
                     eq=eq,
                     ge=ge,
                     goal=goal)

    set_optimizer(model, solver)

    if !verbose
        set_silent(model)
    end

    optimize!(model)

    return objective_value(model)
end



npa_max(expr, level; kw...) = npa_opt(expr, level; goal=:maximise, kw...)
npa_min(expr, level; kw...) = npa_opt(expr, level; goal=:minimise, kw...)
