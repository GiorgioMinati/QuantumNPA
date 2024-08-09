# Assigning party zero to commutative operators is a quick way to ensure commutation
# with all other monomials with minimal change to the rest of the library.

# Commutative operators are identified only by their name
# and ordered lexicographically

abstract type CommutativeOperator <: Operator end

struct Commutative <: CommutativeOperator
    name::String
end

Base.hash(c::CommutativeOperator, h::UInt) = hash(c.name, h)
Base.:(==)(d::CommutativeOperator, c::CommutativeOperator) = d.name == c.name
Base.isless(d::CommutativeOperator, c::CommutativeOperator) = d.name < c.name
Base.string(c::CommutativeOperator) = c.name

function Base.string(c::CommutativeOperator, party::PartyVec) 
    if party != [0]
        @error "Invalid party $party for scalar variable."
    end
    return c.name
end

function is_commutative(optuple::Tuple{PartyVec,Vector{Operator}})
    (p, ops) = optuple
    if p != [0]
        return false
    end
    for o in ops
        if !(typeof(o) <: CommutativeOperator)
            return false
        end
    end
    return true 
end
 
function is_commutative(opv::OpVector)
    for o in opv
        if !is_commutative(o)
            return false
        end
    end
    return true
end

function is_commutative(m::Monomial)
    return is_commutative(m.word)
end

function Base.:*(c::CommutativeOperator, d::CommutativeOperator)
    return c < d ? (1, [c, d]) : (1, [d, c])
end

# For commutative operators we can just them by simply ordering them lexicographically.
function join_commutative(opsx::Vector{Operator}, opsy::Vector{Operator})
   return (1, sort(vcat(opsx, opsy))) 
end

# Commutative monomial constructors
function Monomial(operator::T) where T <: CommutativeOperator
    return Monomial(0, operator)
end

function commutative(name::String)
    return Monomial(Commutative(name))
end

# TODO: add more constructors
# - From list of strings
# - a macro

# TODO: More types of commutative operators
# a macro to create them would be nice
