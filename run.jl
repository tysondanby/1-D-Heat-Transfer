using ForwardDiff, LinearAlgebra, Plots, Interpolations
include("structs.jl")
include("basefunctions.jl")
include("onetoonefunctions.jl")
include("meshing.jl")
include("solving.jl")
include("visualization.jl")

function k(x) 
    if x < 0.03
        return 15
    else
        return 137*exp(25*x-2)
    end
end

ɛ = 0.8
h = 1000
Tinf = 300
Tsur = 250.0
n = 21
Tinit = Tinf
basicmesh = meshingsettings('B',n,linear)
L = .05
Area(x) = 1.0
convectionsourcefunc(T) = h*(Tinf-T)
radiationsourcefunc(T) = -ɛ*5.67E-8*(T^4-Tsur^4)
generationsourcefunc(T) = 4e6
sources = [source((0,0.03),generationsourcefunc)]#[source((0.25,.75),convectionsourcefunc),source((0,1),radiationsourcefunc)]
BCs = [flux(0),convective(h,Tinf)]
layers = [layer("A",0.03),layer("B",0.02)]
testscene=oneDscene(basicmesh,L,Area,sources,BCs,layers,250.0,k)

meshedscene = mesh(testscene)

itterate_solve!(meshedscene) #needed in case source terms are nonlinear
p1=basic_temp_plot(meshedscene)
#p2,p3 = convergencestudy(testscene) #TODO, not working