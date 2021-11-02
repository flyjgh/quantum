module Custom

export f_all, f_line, bellstate, invbellstate, stairs

using BS
using CuYao
using YaoExtensions
using YaoPlots


#------------------------------------------------------------------------

"""
    (n,ƒ) -> chain(n, (k => ƒ) |> put for k ∈ 1:n) |> cu

Put one instance of `ƒ` on every qbit
"""
f_all = (n,ƒ) -> chain(n, (k => ƒ) |> put for k ∈ 1:n) |> cu

"""
    (l,ƒ) -> put(l => ƒ)
    (lines::Union{Vector,Tuple}, ƒ) = put(line => ƒ) for line in lines

`l::Union{Integer, Vector, Tuple, AbstractRange}`

Put `ƒ` on the selected line(s)
"""
f_line(l::Integer,ƒ) = put(l => ƒ)
f_line(lines::Union{Vector,Tuple,AbstractRange},ƒ) = repeat(ƒ, lines)


"""
    bellstate(nb_of_qbits, pair::Pair)
    put(pair.first => H), control(pair.first, pair.second => X)
"""
bellstate(p::Pair) = chain(put(p[1]=>H), control(p[1], p[2]=>X))

"""
    bellstate(nb_of_qbits, pair::Pair)
    control(pair.first, pair.second => X), put(pair.first => H)
"""
invbellstate(p::Pair) = chain(control(p[1], p[2]=>X), put(p[1]=>H))

tp = (n,k,ƒ) -> chain(n, (j=>ƒ) |> put for j ∈ k:n)
stairs = (n,ƒ) -> chain(tp(n,k,ƒ) for k ∈ 1:n) |> cu

end
