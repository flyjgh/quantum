
using BS
import CUDA: CuArray, allowscalar, @allowscalar
import Base.*

allowscalar(false)

mutable struct CuQbit
    state::CuArray{ComplexF64}
    CuQbit() = new(φ0)
    CuQbit(φ::Vector{Number}) = new(φ+0im |> cuda)
    CuQbit(φ::CuArray{ComplexF64}) = new(φ)
end

cuda(q::CuQbit) = CuQbit(CuArray(q.state))
cuda(v::AbstractArray) = CuArray(v)
φ0, φ1 = [1.,0.] .+ 0im |> cuda, [0.,1.] .+ 0im |> cuda

p(q::CuQbit) = @allowscalar [(q.state'φ0)^2,(q.state'φ1)^2]
Base.kron(q::CuQbit,m::AbstractMatrix) = @allowscalar kron(q.state,m)
Base.kron(q::CuQbit,q1::CuQbit) = @allowscalar(q.state ⊗ q1.state)
function *(q::CuQbit,m::AbstractMatrix) @allowscalar(q.state = m*q.state); return q end
*(v::Vector{CuQbit},m::AbstractMatrix) = @allowscalar(m * (v[1] ⊗ v[2]))

σₓ = [0 1
      1 0] .+ 0im |> cuda

σᵧ = im*[0 -1
         1  0] |> cuda

σ𝑧 = [1  0
      0 -1] .+ 0im |> cuda

hadamard = (1 / √2) * [1  1
                       1 -1] .+ 0im |> cuda

cnot = [1 0 0 0
        0 1 0 0
        0 0 0 1
        0 0 1 0] .+ 0im |> cuda

# q = CuQbit()
# q1 = CuQbit()
# q1 * σₓ
# q * hadamard * hadamard * hadamard * hadamard
# [q, q1] * cnot
