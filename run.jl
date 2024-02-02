using ForwardDiff, LinearAlgebra, Plots
include("structs.jl")
include("basefunctions.jl")
include("onetoonefunctions.jl")
include("meshing.jl")
include("solving.jl")

k(x) = 50.0

ɛ = 0.8
h = 10
Tinf = 650
Tsur = 250.0
n = 70
Tinit = Tinf
basicmesh = meshingsettings('B',n,linear)
L = 1.0
Area(x) = 1.0
convectionsourcefunc(T) = h*(Tinf-T)
radiationsourcefunc(T) = -ɛ*5.67E-8*(T^4-Tsur^4)
sources = [source((0.25,.75),convectionsourcefunc),source((0,1),radiationsourcefunc)]
BCs = [flux(0),flux(0)]
layers = [layer("copper",L)]
testscene=oneDscene(basicmesh,L,Area,sources,BCs,layers,250.0,k)

meshedscene = mesh(testscene)

set_up!(meshedscene)

Amatrix, bcolumn = getA_b(meshedscene)
itterate_solve!(meshedscene)
Ts = []
xs = []
for i = 1:1:length(meshedscene.mesh.nodes)
    push!(xs,meshedscene.mesh.nodes[i].pos)
    push!(Ts,meshedscene.mesh.nodes[i].T)
end
plot(xs,Ts)