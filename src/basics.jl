using Yao
using YaoPlots

# |0> = [1. + 0im, 0. + 0im]
# |1> = [0. + 0im, 1. + 0im]
# |+> = (|0> +  |1>) / √2
# |-> = (|0> + i|1>) / √2

circuit = zero_state(1, nbatch=2)
circuit = circuit |> put(1, 1=>H)
circuit.state

1/√2
id =  [1 0
       0 1]

σₓ =  [0 1
      1 0]

σᵧ = [0 -im
     im  0]

σ𝑧 =  [1  0
       0 -1]

S = [1  0
     0 im]

T = √S        #[1 0
              # 0 ℯ^(im*(π/4))]

hadamard = (1 / √2) * [1  1
                       1 -1]

Rₓ = [cos(1/2)    -sin(1/2)im
     -sin(1/2)im   cos(1/2)]

Rᵧ = [cos(1/2) -sin(1/2)
     sin(1/2)  cos(1/2)]

R𝑧 = [ℯ^(1/2im)  0
      0         ℯ^(1/2im)]

cₓ =  [1 0 0 0
      0 1 0 0
      0 0 0 1
      0 0 1 0]

cᵧ =  [1   0   0   0
      0   1   0   0
      0   0   0  -im
      0   0   im  0]

c𝑧 =  [1  0  0  0
       0  1  0  0
       0  0  1  0
       0  0  0 -1]

swap = [1 0 0 0
        0 0 1 0
        0 1 0 0
        0 0 0 1]

ising = [cos(1)    0        0       -im*sin(1)
         0         cos(1)  -sin(1)im 0
         0        -sin(1)im cos(1)   0
        -sin(1)im  0        0        cos(1)]

ccₓ =  [1 0 0 0 0 0 0 0
       0 1 0 0 0 0 0 0
       0 0 1 0 0 0 0 0
       0 0 0 1 0 0 0 0
       0 0 0 0 1 0 0 0
       0 0 0 0 0 1 0 0
       0 0 0 0 0 0 0 1
       0 0 0 0 0 0 1 0]

cswap = [1 0 0 0 0 0 0 0
         0 1 0 0 0 0 0 0
         0 0 1 0 0 0 0 0
         0 0 0 1 0 0 0 0
         0 0 0 0 1 0 0 0
         0 0 0 0 0 0 1 0
         0 0 0 0 0 1 0 0
         0 0 0 0 0 0 0 1]

bellstate = cₓ * (hadamard ⊗ id)

invbell = (hadamard ⊗ id) * cₓ

rotationₓ(θ) = [cos(θ/2)   -sin(θ/2)im
               -sin(θ/2)im  cos(θ/2)]

rotationᵧ(θ) = [cos(θ/2) -sin(θ/2)
                sin(θ/2)  cos(θ/2)]

rotation𝑧(θ) = [ℯ^(θ/2im) 0
                0        ℯ^(θ/2im)]

# q = CuQbit()
# q1 = CuQbit()
# q1 * σₓ
# q * hadamard * hadamard * hadamard * hadamard
# [q, q1] * cnot
