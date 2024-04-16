using ForwardDiff, LinearAlgebra, Plots, Interpolations
include("structs.jl")
include("basefunctions.jl")
include("onetoonefunctions.jl")
include("meshing.jl")
include("solving.jl")
include("visualization.jl")

#=
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
=#
#=
function homework5()
    function Γ(x)
        return 1
    end
    function get_velocity_func(v)
        function  func(x)
            return [v,0,0]
        end
        return func
    end
    function rho(x)
        return 1
    end
    function A(x)
        return 1
    end
    n = 21
    L = 1
    Tsur = 1
    mesh_settings = meshingsettings('A',n,linear)
    BCs = [constanttemp([0.0,0.0,0.0],1),constanttemp([1.0,0.0,0.0],0)]
    layers = [layer("A",L)]
    scenes = []
    for velocity in [1.0,20.0,75.0]
        for scheme in ["central","upwind","hybrid","power"]
            push!(scenes,oneDscene(mesh_settings,scheme,L,A,[],BCs,layers,Tsur,Γ,get_velocity_func(velocity),rho))
        end
    end
    meshedscenes = []
    for scene in scenes
        push!(meshedscenes,mesh(scene))
    end
    
    for scene in meshedscenes
        itterate_solve!(scene)
    end

    println()
    println("PROBLEM 1 ---------------------------------")
    println()
    velocities =[1.0,20.0,75.0]
    plots = []
    tableseries = []
    push!(tableseries,["Central","Upwind","Hybrid","Power Law"])
    for k = 1:1:length(velocities)
        profiles = []
        velocity = velocities[k]
        labels = ["Central","Upwind","Hybrid","Power"]
        xs = []
        exactprofile = []
        series_of_errors = []
        for j = 1:1:length(labels)
            index = (k-1)*4+j
            s = meshedscenes[index]
            xs = []
            ϕs = []
            error = 0.0
            for i = 1:1:length(s.mesh.nodes)
                push!(xs,s.mesh.nodes[i].pos[1])
                if j == 1
                    push!(exactprofile,1 - (exp(xs[i]*velocity)-1)/(exp(velocity)-1))
                end
                push!(ϕs,s.mesh.nodes[i].T)
                error = error + abs(ϕs[i]-exactprofile[i])
            end
            push!(series_of_errors,error)
            push!(profiles,ϕs)
        end
        push!(tableseries,series_of_errors)
        push!(profiles,exactprofile)
        println()
        println()
        println("Table of values of ϕ for all schemes and u = $velocity")
        println()
        printtable([xs,profiles...],["x",labels...,"Exact"],decimals = 4)

        push!(plots,plot(xs,profiles,xlabel = "Distance", ylabel = "ϕ",labels=reshape([labels...,"Exact"],(1,5)),title = "ϕ for All Schemes and u = $velocity"))
    end
    println()
    println()
    println("Table of The Error of Each Approcimate Solution")
    println()
    printtable(tableseries,["Scheme","u = 1","u = 20","u = 75"],decimals = 4)
    return plots
end
=#
function homework6()
    Tsur = 293.0
    Tb = 353.0
    Tinf = 293.0
    σ = 5.67E-8
    w = 1.0
    b = 0.0003
    n = 21
    L = .01
    function k(x)
        return 177
    end
    function get_velocity_func(v)
        function  func(x)
            return [v,0,0]
        end
        return func
    end
    function rho(x)
        return 2770.0
    end
    function A(x)
        return w*b
    end
    function Tinit(pos)
        return 293.0
    end
    function cp(x)
        return 875.0
    end
    function knownT(xin)
        x = xin[1]
        h = 30.0
        P = 2*w+2*b
        m = sqrt(h*P/(k(x)*A(x)))
        theta = cosh(m*(L-x))/cosh(m*L)
        return theta*(Tb-Tinf) + Tinf
    end

    plots = []
    #--------Perform grid independence study on emissivity 0 and h 30
    h = 30.0
    ɛ = 0.0
    function S1(T)
        return -(2/(b))*(h*(T-Tinf)+ɛ*σ*(T^4-Tinf^4))
    end
    mesh_settings = meshingsettings('A',n,linear)
    BCs = [constanttemp([0.0,0.0,0.0],Tb),flux([1.0,0.0,0.0],0)]
    layers = [layer("A",L)]
    testscene = oneDscene(mesh_settings,"central",L,A,[source((0,L),S1)],BCs,layers,Tsur,k,get_velocity_func(0.0),rho,cp,Tinit)
    push!(plots,convergencestudy_grid(testscene,(10,400)))#convergencestudy_steady(testscene))
    meshedtestscene = mesh(testscene)
    itterate_solve_steady!(meshedtestscene)
    global p1 = basic_temp_plot(meshedtestscene)
    println("2A done")
    n = 300#MANUALLY INPUT THE NUMBER OF NODES NEEDED FOR CONVERGENCE

    dx = (L)/(n-1)
    xsteps = collect(range(0,L,n))
    mesh_settings = meshingsettings('A',n,linear)
    scenes = []
    Qhistories = []
    for h in [5,30]
        for ɛ in [0.0,1.0]
            function S(T)
                return -(2/b)*(h*(T-Tinf)+ɛ*σ*(T^4-Tsur^4))
            end
            push!(scenes,oneDscene(mesh_settings,"central",L,A,[source((0,L),S)],BCs,layers,Tsur,k,get_velocity_func(0.0),rho,cp,Tinit))
        end
    end
    meshedscenes = []
    for scene in scenes
        push!(meshedscenes,mesh(scene))
    end
    trange = (0,5)
    #push!(plots,convergence_study_time(deepcopy(scenes[3]),trange,(100,3000))) 
    println("2B done")

    dt = 0.0025 #MANUALLY INPUT THE dt NEEDED FOR CONVERGENCE


    timesteps = collect(trange[1]:dt:trange[2])
    for i = 1:1:length(meshedscenes)
        scene = meshedscenes[i]
        Thistory = unsteady_solve!(scene,timesteps)
        Qs = []
        for Ts in Thistory
            dTdx = (Ts[2]-Ts[1])/dx
            push!(Qs,-k(0.0)*w*b*dTdx)
        end
        h = 0.0
        ɛ = 0.0
        if i == 1
            h=5.0
            ɛ = 0.0
        elseif i == 2
            h = 5
            ɛ = 1.0
        elseif i ==3
            h = 30.0
            ɛ = 0.0
        else
            h = 30.0
            ɛ = 1.0
        end
        function knownTs(xin)
            x = xin[1]
            P = 2*w+2*b
            m = sqrt(h*P/(k(x)*w*b))
            theta = cosh(m*(L-x))/cosh(m*L)
            return theta*(Tb-Tinf) + Tinf
        end
        
        push!(plots,plot(timesteps,Qs,xlabel = "Time (s)", ylabel = "Q(t) at Base of Fin (W)",title = "Heat Transfer Over Time for h = $h and ɛ = $ɛ"))
        push!(plots,plot(xsteps,[Thistory[end],@. knownTs(xsteps)],xlabel = "Position (m)", ylabel = "Temperature (K)",labels=reshape(["Numerical Method","Exact"],(1,2)),title = "St-St Temp Profiles for h = $h and ɛ = $ɛ"))
    end

    return plots
end

function homework4()
    ɛ = 1.0
    h = 10
    Tinf = 273
    Tb = 400
    Ts = 273
    n = 82
    basicmeshes = []
    function k(x) 
        return 401
    end
    function fT(x)
        return 273.15
    end
    function Teven(x)
        return 400.0
    end
    function getzero(x)
        return [0.0, 0.0, 0.0]
    end
    L = .02
    D = .003

    Area(x) = pi*D^2
    sourcefunc(T) = -(4/D)*(h*(T-Tinf)+ɛ*5.67E-8*(T^4-Ts^4))#/Area(0.0)
    sources = [source((0,L),sourcefunc)]#[source((0.25,.75),convectionsourcefunc),source((0,1),radiationsourcefunc)]
    BCs = [constanttemp([0.0,0.0,0.0],Tb),flux([L,0.0,0.0],0.0)]
    layers = [layer("A",L)]

    en = sqrt(4*h/(D*k(0)))
    solution(x) = (cosh(en*(L-x))/cosh(en*L))*(Tb-Tinf) + Tinf
    basicmesh = meshingsettings('B',n,linear)
    scene=oneDscene(basicmesh,"central",L,Area,sources,BCs,layers,Tinf,k,getzero,Teven,Teven,fT)
    meshedscene = mesh(scene)
    xs = []
    startTs = []
    for node in meshedscene.mesh.nodes
        push!(xs,node.pos[1])
        push!(startTs,solution(node.pos[1]))
    end
    #println(xs)
    println("no errors until solving")
    itterate_solve_GS!(meshedscene,startTs)
    global p3=basic_temp_plot(meshedscene)
    exs = collect(0:.001:L)
    global sln = @. solution(exs)
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