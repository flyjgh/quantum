# using Yao
# using YaoPlots
# using Plots

# |0> = [1. + 0im, 0. + 0im]
# |1> = [0. + 0im, 1. + 0im]
# |+> = (|0> + |1>) / √2 == [1/√2 + 0im, 1/√2 + 0im]
# |_> = (|0> - |1>) / √2 == [1/√2 + 0im,-1/√2 + 0im]

# X * |0> == |1>
# H * |0> == |+>
# H * |1> == |_>

# A state |ψ> == α|0> + β|1> == 1
# probability == |ψ|²

# Hadamard
# H(α|0> + β|1>) = α/√2(α|0> + β|1>) + β/√2(α|0> - β|1>)
#                = (α+β)/√2  * (|0>) + (α-β)/√2  * (|1>)

# circuit = zero_state(1, nbatch=2)
# circuit = circuit |> put(1, 1=>H)
# circuit.state

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

H = (1/√2) * [1  1
              1 -1] .+ 0im

RX = [cos(1/2)    -sin(1/2)im
     -sin(1/2)im   cos(1/2)]

RY = [cos(1/2) -sin(1/2)
      sin(1/2)  cos(1/2)]

RZ = [ℯ^(1/2im)  0
      0          ℯ^(1/2im)]

cx =  [1 0 0 0
       0 1 0 0
       0 0 0 1
       0 0 1 0]

cy =  [1   0   0   0
       0   1   0   0
       0   0   0  -im
       0   0   im  0]

cz =  [1  0  0  0
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

ccX =  [1 0 0 0 0 0 0 0
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

bell = cx * (H ⊗ id)

invbell = (H ⊗ id) * cx

rotationX(θ) = [cos(θ/2)   -sin(θ/2)im
               -sin(θ/2)im  cos(θ/2)]
rotationX(x, θ) = rotationX(θ) * x

rotationY(θ) = [cos(θ/2) -sin(θ/2)
                sin(θ/2)  cos(θ/2)]
rotationY(x, θ) = rotationY(θ) * x

rotationZ(θ) = [ℯ^(θ/2im) 0
                0        ℯ^(θ/2im)]
rotationZ(x, θ) = rotationZ(θ) * x

phase_shift(θ) = [1 0
                 0 ℯ^(θ*im)]
phase_shift(x, θ) = phase_shift(θ) * x

φ0 = [1. + 0im, 0. + 0im]
φ1 = [0. + 0im, 1. + 0im]
φx = [1/√2 + 0im, 1/√2 + 0im] # == H * φ0
φy = [1/√2 + 0im,-1/√2 + 0im] # == H * φ1

σₓ * φ0 == φ1
H * φ0 == φx
H * φ1 == (1/√2)*(φ0-φ1) == φy
φy == (φ0 - φ1) / √2 == [1/√2 + 0im,-1/√2 + 0im]

H * φ0
swap * (φ1 ⊗ φ0)

invbell * bell * (φ1 ⊗ φ0) ≈ (φ1 ⊗ φ0)
