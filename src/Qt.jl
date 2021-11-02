
#-----------------------------------------------------------------------

module Qt

export Qbit, ⊗, φ0, φ1, φx, φy, randstate, apply, proba, kron, merge, *
export id, σₓ, σᵧ, σ𝑧, Rₓ, Rᵧ, R𝑧, cₓ, cᵧ, c𝑧, ccₓ, S, T, swap, cswap
export hadamard, ising, bellstate, invbell
export rotationₓ, rotationᵧ, rotation𝑧, states

#-----------------------------------------------------------------------

import Base.*

const φ0 = [1.,0.] .+ 0im
const φ1 = [0.,1.] .+ 0im
const φx = (φ0 + im*φ1) / √2
const φy = (φ0 + im*φ1) / √2
const ⊗ = kron

randstate(nb_qbit) = im .* rand(2^nb_qbit)

mutable struct Qbit
    state::Vector{ComplexF64}
    n::Int64
    Qbit() = new(φ0,1)
    Qbit(args::Vararg{Any,N}) where {N} = [Qbit.(args)...]
    Qbit(q::Qbit) = q
    Qbit(v::Vector{Qbit}) = merge(v)
    Qbit(φ::Vector) = new(ComplexF64.(φ), length(φ) |> log2)
end

#-----------------------------------------------------------------------

states(n) = @inbounds [ComplexF64.(BitArray(v==x for v∈1:n)) for x∈1:n]

proba(q::Qbit) = @inbounds [(q.state'state)^2 for state ∈ states(length(q.state))]

apply(v::Vector{Qbit}, m::Matrix{T}) where T = vcat([q for q ∈ v]) .* Ref(m)

Base.kron(q::Qbit,m::Matrix{T}) where T = kron(q.state,m)

Base.kron(q::Qbit,q1::Qbit) = q.state ⊗ q1.state

Base.kron(v::Vector{Qbit}) = @inbounds ⊗([q.state for q in v]...)

Base.merge(v::Vector{Qbit}) = Qbit(⊗(v))

*(v::Vector{Qbit}, m::Matrix{T}) where T = Qbit.(m * ⊗(v))

function *(q::Qbit,m::Matrix{T}) where T  q.state = m * q.state ; q end

#-----------------------------------------------------------------------

const id =  [1 0
             0 1]

const σₓ =  [0 1
             1 0]

const σᵧ=im*[0 -1
             1  0]

const σ𝑧 =  [1  0
             0 -1]

const S = [1  0
           0 im]

const T = √S        #[1 0
                    # 0 ℯ^(im*(π/4))]
const hadamard = (1 / √2) * [1  1
                             1 -1]

const Rₓ = [cos(1/2) -sin(1/2)im
           -sin(1/2)im  cos(1/2)]

const Rᵧ = [cos(1/2) -sin(1/2)
            sin(1/2)  cos(1/2)]

const R𝑧 = [ℯ^(1/2im) 0
            0 ℯ^(1/2im)]

const cₓ =  [1 0 0 0
             0 1 0 0
             0 0 0 1
             0 0 1 0]

const cᵧ =  [1 0 0   0
             0 1 0   0
             0 0 0 -im
             0 0 im  0]

const c𝑧 =  [1 0 0  0
             0 1 0  0
             0 0 1  0
             0 0 0 -1]

const swap = [1 0 0 0
              0 0 1 0
              0 1 0 0
              0 0 0 1]

const ising = [cos(1) 0 0 -im*sin(1)
              0 cos(1)  -sin(1)im 0
              0 -sin(1)im  cos(1) 0
              -sin(1)im  0 0 cos(1)]

const ccₓ =  [1 0 0 0 0 0 0 0
              0 1 0 0 0 0 0 0
              0 0 1 0 0 0 0 0
              0 0 0 1 0 0 0 0
              0 0 0 0 1 0 0 0
              0 0 0 0 0 1 0 0
              0 0 0 0 0 0 0 1
              0 0 0 0 0 0 1 0]

const cswap = [1 0 0 0 0 0 0 0
               0 1 0 0 0 0 0 0
               0 0 1 0 0 0 0 0
               0 0 0 1 0 0 0 0
               0 0 0 0 1 0 0 0
               0 0 0 0 0 0 1 0
               0 0 0 0 0 1 0 0
               0 0 0 0 0 0 0 1]

const bellstate = cₓ * (hadamard ⊗ id)

const invbell = (hadamard ⊗ id) * cₓ

rotationₓ(θ) = [cos(θ/2) -sin(θ/2)im
               -sin(θ/2)im  cos(θ/2)]

rotationᵧ(θ) = [cos(θ/2) -sin(θ/2)
                sin(θ/2)  cos(θ/2)]

rotation𝑧(θ) = [ℯ^(θ/2im) 0
                0 ℯ^(θ/2im)]

#-----------------------------------------------------------------------

end
