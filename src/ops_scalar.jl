# WIP: Scalar extension support
#  

# A scalar factor is just a commutative variable with an associated monomial.

struct ScalarFactor <: CommutativeOperator
    name::String
    monomial::Monomial
end

# Scalar factor constructors 
function ScalarFactor(m::Monomial)
    return ScalarFactor("<$m>", m)
end

function scalarfactor(name::String, m::Monomial)
    return Monomial(ScalarFactor(name, m))
end

function scalarfactor(m::Monomial)
    return Monomial(ScalarFactor(m))
end
# TODO: add some macro

# Check if a monomial contains a single scalar factor
function is_scalar_factor(m::Monomial)
    w = m.word
    return length(w) == 1 && is_commutative(w) && 
           length(w[1][2]) == 1 && typeof(w[1][2][1]) == ScalarFactor  
end

# Connectivity graph that determines when we can substitute the scalar factors.
# It is a list edges on the party numbers.
ConnectionStructure = Set{Tuple{Int, Int}}

Base.sort(t::Tuple{Int, Int}) = t[1] > t[2] ? (t[2], t[1]) : t

function connection_structure(edges::Tuple{Int, Int}...)
    edges = map(sort, edges)
    edges = filter(t->t[1]!=t[2], edges)
    return ConnectionStructure(edges)       
end

function connection_structure(con_parties::Tuple{String, String}...)
    con_nums = map(x::Tuple{String, String}->map(party_num, x), con_parties)
    return connection_structure(con_nums...)
end

# It is often easier to define the complementary version of this graph
function separability_structure(sep_edges::Tuple{Int, Int}...)
    con = ConnectionStructure()
    all_parties = collect(1:maximum([(sep_edges...)...]))
    unique!(all_parties)
    n = length(all_parties)
    for i=1:n, j=i+1:n
        p1, p2 = all_parties[i], all_parties[j]
        if p1 <= 0 || p2 <= 0
            @error "Invalid party identifier"
        end  
        if !((p1, p2) in sep_edges || (p2, p1) in sep_edges)
            push!(con, (p1, p2))
        end
    end
    return con
end

function separability_structure(sep_parties::Tuple{String, String}...)
    sep_nums = map(x::Tuple{String, String}->map(party_num, x), sep_parties)
    return separability_structure(sep_nums...)
end

function is_sep(pvec1::PartyVec, pvec2::PartyVec, con::ConnectionStructure)
    if pvec1 == [0] || pvec2 == [0]
        return true
    end
    for p1 in pvec1, p2 in pvec2
        if (p1, p2) in con || (p2, p1) in con
            return false
        end
    end
    return true
end

function _extract_connected!(ops::OpVector, con::ConnectionStructure)
    connected_ops = [popfirst!(ops)]
    connected_pvec = [connected_ops[1][1]...]
    all_sep = false
    while !all_sep && length(ops) > 0
        all_sep = true
        for (i, op) in enumerate(ops)
            if !is_sep(connected_pvec, op[1], con)
                push!(connected_pvec, op[1]...)
                pushfirst!(connected_ops, op)
                remin!(connected_ops)
                deleteat!(ops, i)
                all_sep = false
            end
        end
    end
    return connected_ops
end

# Split monomials according to connection rules
function split_monomials(m::Monomial, con::ConnectionStructure)
    split_mons = Vector{Monomial}()
    ops = copy(m.word)
    while length(ops) > 0
        connected_ops = _extract_connected!(ops, con)
        push!(split_mons, Monomial(connected_ops))
    end
    return split_mons
end

function split_reduce_expr(mon::Monomial, con::ConnectionStructure, space::Linspace)
    smon = split_monomials(mon, con)
    new = Id
    for m in smon
        m = reduce_expr!(m, space)
        new *= m
    end
    return new
end

function split_reduce_expr(poly::Polynomial, con, space::Linspace)
    new_poly = Polynomial(poly.cfsize, poly.blockstruct)
    for (mon, coeff) in collect(poly.terms)
        new_mon = split_reduce_expr(mon, con, space)
        new_poly += new_mon*coeff
    end
    return new_poly
end

# ScalarExtension
# 
# This struct contains all the information about the scalar extension
# and the associated methods allow to build it quickly.
# The factors are stored in a dictionary indexed by monomials for 
# quick access (probably not needed).

ScalarFactors = Dict{Monomial,Monomial}

struct ScalarExtension
    factors::ScalarFactors
    graph::ConnectionStructure
end

function Base.getindex(se::ScalarExtension, m::Monomial)
    if haskey(se.factors, m)
        return se.factors[m]
    else
        return m
    end
end

# Functions to add scalar factors
function push_factor!(se::ScalarExtension, f::ScalarFactor)
    if !haskey(se.factors, f.monomial)
        se.factors[f.monomial] = Monomial(f)
    end 
end

function push_factor!(se::ScalarExtension, m::Monomial)
    if is_scalar_factor(m)
        f = m.word[1][2][1]
        se.factors[f.monomial] = Monomial(f)
    end
end

function push_factors!(se::ScalarExtension, items...)
    for item in items
        push_factor!(se, item)
    end
end

# Functions to add scalar factors created from monomial
function new_factor!(se::ScalarExtension, m::Monomial)
    mf = ScalarFactor(m)
    push_factor!(se, mf)
end

function new_factors!(se::ScalarExtension, mons...)
    for m in mons
        new_factor!(se, m)
    end
end

# Functions to add edges
function add_edge!(se::ScalarExtension, e::Tuple{Int, Int})
    if e[1] == e[2]
        @error "Loop edges are not allowed."
    end
    e = sort(e)
    push!(se.graph, e)
end

function add_edges!(se::ScalarExtension, edges...)
    for e in edges
        add_edge!(se, e)
    end
end

# Constructors
function ScalarExtension(factors, con::ConnectionStructure)
    se = ScalarExtension(ScalarFactors(), con)
    push_factors!(se, factors...)
    return se
end

function ScalarExtension(factors, edges, connection_edges=false)
    if connection_edges
        con = connection_structure(edges...)
    else
        con = separability_structure(edges...)
    end
    se = ScalarExtension(ScalarFactors(), con)
    push_factors!(se, factors...)
    return se
end

function define_scalar_extension(mons, con::ConnectionStructure)
    se = ScalarExtension(ScalarFactors(), con)
    new_factors!(se, mons...)
    return se
end

function define_scalar_extension(mons, edges, connection_edges=false)
    if connection_edges
        con = connection_structure(edges...)
    else
        con = separability_structure(edges...)
    end
    se = ScalarExtension(ScalarFactors(), con)
    new_factors!(se, mons...)
    return se
end

# Functions to perform scalar factors substitutions when appropriate
#

function subs_scalar_factors(mon::Monomial, se::ScalarExtension; eq::Linspace=Linspace())
    smon = split_monomials(mon, se.graph)
    new_mon = Id
    for m in smon
        if m in keys(se.factors)
            new_mon *= se[m]
        else
            new_mon *= m
        end
    end
    return new_mon
end

function subs_scalar_factors(poly::Polynomial, se::ScalarExtension; eq::Linspace=Linspace())
    new_poly = Polynomial(poly.cfsize, poly.blockstruct)
    for (mon, coeff) in collect(poly.terms)
        new_mon = subs_scalar_factors(mon, se, eq=eq)
        new_poly += new_mon*coeff
    end
    return new_poly
end

# Reduce polynomials taking into account scalar extension factors

_lappend_monomial(space::Linspace, m::Monomial) = Linspace(m*k => m*v for (k,v) in space)

function reduce_se_expr!(expr, space::Linspace, se::ScalarExtension)
    factors = collect(values(se.factors))
    for f in factors
        fspace = _lappend_monomial(space, f)
        reduce_expr!(expr, fspace)
        ffspace = _lappend_monomial(fspace, f)
        reduce_expr!(expr, ffspace)
    end    
end

function reduce_se_expr(expr, space::Linspace, se::ScalarExtension)
    new_expr = new_polynomial(expr)    
    reduce_se_expr!(new_expr, space, se)
    return new_expr
end

# Overload the NPA builder functions to construct the moment matrix with scalar factors 
#

function npa_moment(operators::Vector{<:Union{Monomial,Polynomial}}, se::ScalarExtension)
    factors = collect(values(se.factors))
    return npa_moment([operators; factors])
end

npa_moment(source, level, se::ScalarExtension) = npa_moment(npa_level(source, level), se::ScalarExtension)

function _reduce_problem(expr, moment::Polynomial, se::ScalarExtension; eq=[], ge=[])
    # Reduce constraints to canonical form
    expr = Polynomial(conj_min(expr))
    eq = linspace(map(conj_min, eq))
    ge = map(x->Polynomial(conj_min(x)), ge)

    if haskey(eq, Id)
        @error "Contradiction Id = 0 in equality constraints."
    end

    # Reduce the objective expression, using constraints to eliminate
    # monomials.
    expr = subs_scalar_factors(expr, se, eq=eq)
    expr = split_reduce_expr(expr, se.graph, eq)

    # Reduce moments using equality constraints.
    moment = subs_scalar_factors(moment, se, eq=eq)
    moment = split_reduce_expr(moment, se.graph, eq)

    # Reduce inequality constraints then include them as inequalities along
    # with the original moment matrix.
    ge = map(x->subs_scalar_factors(x, se, eq=eq), ge)
    ge = reduce_exprs(ge, eq)
    
    return expr, moment, ge, eq
end

# Define SDP for scalar variables
function npa2sdp(expr, moment::Polynomial, se::ScalarExtension; eq=[], ge=[])
    expr, moment, ge, eqspace = _reduce_problem(expr, moment, se, eq=eq, ge=ge) 
    

    #reduce_se_expr!(expr, eqspace, se)
    #reduce_se_expr!(moment, eqspace, se)
    #ge = map(x->reduce_se_expr(x, eqspace, se), ge)

    return (expr, vcat([moment], ge))
end

function npa2sdp(expr, level, se::ScalarExtension; eq=[], ge=[])
    moment = npa_moment([expr, eq, ge], level, se)
    return npa2sdp(expr, moment, se, eq=eq, ge=ge)
end

function npa2jump(expr, level_or_moments, se::ScalarExtension;
                  eq=[],
                  ge=[],
                  goal=:maximise,
                  solver=nothing,
                  verbose=nothing)
    (expr, moments) = npa2sdp(expr, level_or_moments, se, eq=eq, ge=ge)

    model = sdp2jump(expr, moments,
                     goal=goal,
                     solver=solver,
                     verbose=verbose)
    return model
end

# TODO: test for intermediate levels like 1+AB
