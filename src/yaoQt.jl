using CuYao
using YaoExtensions
using YaoPlots
include("Customgates.jl")
using .Custom

# begin
# 	bellcircuit = chain(2, put(1=>H), control(1, 2=>X))
# 	plot(bellcircuit)
# end
# q1 = ArrayReg(bit"00")
# a= q1 |> bellcircuit
# a |> state

invbell = invbellstate(2 => 1)
bell = bellstate(1 => 2)
h  = f_all(5,H)
x  = put(4 => X)
xy = f_line([2,5],chain(H,X))
cx = control(1, 2=>X)
hz = f_line(1:3,Y)
c=chain(2,bell,invbell)
c |> plot
c.blocks
tmpcircuit = zero_state(2) |> cu
circuit = tmpcircuit |> c
circuit.state

A(i, j) = control(i, j=>shift(2π/(1<<(i-j+1))))
B(n, k) = chain(n, j==k ? put(k=>H) : A(j, k) for j in k:n)
qft(n) = chain(B(n, k) for k in 1:n) |> cu

cct = qft(3) |> cu
plot(cct)
cct |> ♯
cct[1]
