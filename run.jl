using ForwardDiff, LinearAlgebra, Plots, Interpolations
include("structs.jl")
include("basefunctions.jl")
include("onetoonefunctions.jl")
include("meshing.jl")
include("solving.jl")
include("visualization.jl")

function homework3_problem1()
    ɛ = 0.8
    h = 1000
    Tinf = 300
    nCVs = [5,10,20,40,80,160,320,640,1280,2560,5120,10240]
    ns = @. nCVs + 2# Number CVs +2
    basicmeshes = []
    function k(x) 
        if x < 0.03
            return 15
        else
            return 137*exp(25*x-2)
        end
    end
    L = .05
    Area(x) = 1.0
    generationsourcefunc(T) = 4e6
    sources = [source((0,0.03),generationsourcefunc)]#[source((0.25,.75),convectionsourcefunc),source((0,1),radiationsourcefunc)]
    BCs = [flux(0),convective(h,Tinf)]
    layers = [layer("A",0.03),layer("B",0.02)]
    scenes = []
    for i = 1:1:length(ns)
        push!(basicmeshes,meshingsettings('B',ns[i],linear))
        push!(scenes,oneDscene(basicmeshes[i],L,Area,sources,BCs,layers,250.0,k))
    end
    
    meshedscenes = @. mesh(scenes)
    @. itterate_solve!(meshedscenes)
    fits = []
    for j = 1:1:length(ns)
        newxs = []
        newTs = []
        for i = 1:1:length(meshedscenes[j].mesh.nodes)
            push!(newTs,meshedscenes[j].mesh.nodes[i].T )
            push!(newxs,meshedscenes[j].mesh.nodes[i].pos)
        end
        push!(fits,linear_interpolation(newxs, newTs))
    end
    location = 0.01
    println()
    println()
    println()
    println("Part a:")
    println("Temperature at location x=$location"*" instead of the left boundary due to the left boundary always giving the exact values at the boundaries for this problem with this code. This effect was thouroghly discussed on my last homework.")
    T1 = fits[1](location)
    println("For 5 CVs: T(x=$location"*") = $T1")
    T2 = fits[2](location)
    println("For 10 CVs: T(x=$location"*") = $T2")
    T3 = fits[3](location)
    println("For 20 CVs: T(x=$location"*") = $T3")
    order, T = richardson_extrapolate_results((T1,T2,T3),(nCVs[1],nCVs[2],nCVs[3]))
    println()
    println("Part b:")
    println("For grid independence: T(x=$location"*") = $T")
    println("Order: $order")
    println()
    println("Part c:")
    println("PLOT 1 HERE")
    
    temps = []
    percentdiffs=[]
    for i = 1:1:length(ns)
        push!(temps,fits[i](location))
        push!(percentdiffs,(temps[i]-T)*100/T)
    end
    global p1 = dedupe_and_correlate(nCVs,temps;xlabel="Number of CVs",ylabel="Temperature evaluated at x=$location")
    plot!([nCVs[1],nCVs[end]],[T,T])
    printtable((nCVs,temps,percentdiffs),("Number of CVs","Temperature evaluated at x=$location","Percent Error"))
end

function homework3_problem2()
    ɛ = 0.0
    h = 10
    Tinf = 273
    Tb = 400
    Ts = 273
    n = 20
    basicmeshes = []
    function k(x) 
        return 401
    end
    L = .02
    D = .003

    Area(x) = 1.0
    sourcefunc(T) = -(4/D)*(h*(T-Tinf)+ɛ*5.67E-8*(T^4-Ts^4))
    sources = [source((0,L),sourcefunc)]#[source((0.25,.75),convectionsourcefunc),source((0,1),radiationsourcefunc)]
    BCs = [constanttemp(Tb),flux(0)]
    layers = [layer("A",L)]

    en = sqrt(4*h/(D*k(0)))
    solution(x) = (cosh(en*(L-x))/cosh(en*L))*(Tb-Tinf) + Tinf
    basicmesh = meshingsettings('B',n,linear)
    scene=oneDscene(basicmesh,L,Area,sources,BCs,layers,250.0,k)
    global p2=convergencestudy_knownsolution(scene,solution;reltol =1e-4)
    basicmesh = meshingsettings('B',n,linear)

    scene=oneDscene(basicmesh,L,Area,sources,BCs,layers,250.0,k)
    meshedscene = mesh(scene)
    itterate_solve!(meshedscene) #needed in case source terms are nonlinear
    global p3=basic_temp_plot(meshedscene)
    exs = collect(0:.001:L)
    sln = @. solution(exs)
    dTdx = (meshedscene.mesh.nodes[2].T-meshedscene.mesh.nodes[1].T)/(meshedscene.mesh.nodes[2].pos-meshedscene.mesh.nodes[1].pos)
    qs = -k(0)*dTdx
    println("qs'' = $qs")
    #plot!(exs,sln)
end

function homework3_problem3()
    ɛ = 1.0
    h = 10
    Tinf = 273
    Tb = 400
    Ts = 273
    n = 47
    basicmeshes = []
    function k(x) 
        return 401
    end
    L = .02
    D = .003

    Area(x) = 1.0
    sourcefunc(T) = -(4/D)*(h*(T-Tinf)+ɛ*5.67E-8*(T^4-Ts^4))
    sources = [source((0,L),sourcefunc)]#[source((0.25,.75),convectionsourcefunc),source((0,1),radiationsourcefunc)]
    BCs = [constanttemp(Tb),flux(0)]
    layers = [layer("A",L)]

    en = sqrt(4*h/(D*k(0)))
    #solution(x) = (cosh(en*(L-x))/cosh(en*L))*(Tb-Tinf) + Tinf
    basicmesh = meshingsettings('B',n,linear)
    scene=oneDscene(basicmesh,L,Area,sources,BCs,layers,250.0,k)
    global p4 =convergencestudy(scene;reltol =1e-6)
    basicmesh = meshingsettings('B',n,linear)

    scene=oneDscene(basicmesh,L,Area,sources,BCs,layers,250.0,k)
    meshedscene = mesh(scene)
    itterate_solve!(meshedscene) #needed in case source terms are nonlinear
    global p5=basic_temp_plot(meshedscene)
    #exs = collect(0:.001:L)
    #sln = @. solution(exs)
    #plot!(exs,sln)
end


#=
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
=#