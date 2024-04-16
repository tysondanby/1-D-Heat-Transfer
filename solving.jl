function record_T(s)
    Ts = []
    for i = 1:1:length(s.mesh.nodes)
        push!(Ts,s.mesh.nodes[i].T)
    end
    return Ts
end


function find_a_b_interior(s,n)#n is the node number
    node = s.mesh.nodes[n]
    function totalsource(T)
        S = 0
        if length(s.sources) >= 1
            for i = 1:1:length(s.sources)
                source = s.sources[i]
                if source.range[1] < node.pos < source.range[2]
                    S = S + source.func(T)
                end
            end
        end
        return S
    end
    kP = s.k(node.pos)
    sp = ForwardDiff.derivative(totalsource,node.T)
    sc = totalsource(node.T) - node.T*sp
    an = []
    for i = 1:1:length(node.neighbors)
        neighbornode = s.mesh.nodes[node.neighbors[i]]
        kN = s.k(neighbornode.pos)
        for boundaryindex in node.boundaries #Find the boundary associated with this neighbor
            boundary = s.mesh.boundaries[boundaryindex]
            if (node.pos <= boundary.pos <= neighbornode.pos ) || (node.pos >= boundary.pos >= neighbornode.pos )
                del_n = abs(node.pos-neighbornode.pos)
                del_n_P = abs(boundary.pos - node.pos)#amount of del_n in the P region
                del_n_N = del_n - del_n_P
                kn = del_n / (del_n_P/kP + del_n_N/kN)
                push!(an,kn/del_n)
            end
        end
    end
    boundary1 = s.mesh.boundaries[node.boundaries[1]]
    boundary2 = s.mesh.boundaries[node.boundaries[end]]
    del_x_p = abs(boundary2.pos - boundary1.pos)#TODO: this is not generalized at all, only works for my cases.
    ap = sum(an) - del_x_p * sp
    b = del_x_p*sc
    return ap,an,b
end

function find_a_b_zero_volume(s,n)
    node = s.mesh.nodes[n]
    neighbornode = s.mesh.nodes[node.neighbors[1]] #TODO: for 1D there should only be one, but it will be complicated in higher dimensions.
    ap = 0
    an = [0]
    b = 0
    kn = s.k(neighbornode.pos)
    del_n = abs(node.pos-neighbornode.pos)
    if n == 1
        n=1
    else 
        n = 2
    end
    BC = s.BCs[n]
    if typeof(BC) <:flux
        ap = kn / del_n
        an = [kn / del_n]
        b = BC.value
    elseif typeof(BC) <:convective
        ap = kn / del_n + BC.value
        an = [kn / del_n]
        b = BC.value*BC.Tinf
    elseif typeof(BC) <:constanttemp
        ap = 1.0
        an=[0]
        b = BC.value
    else
        error("Invalid boundary condition type.")
    end
    return ap,an,b
end

function find_a_b(s,n)
    node = s.mesh.nodes[n]
    function totalsource(T)
        S = 0.0
        if length(s.sources) >= 1
            for i = 1:1:length(s.sources)
                source = s.sources[i]
                if source.range[1] < node.pos[1] < source.range[2]
                    S = S + source.func(T)
                end
            end
        end
        return S
    end
    kP = s.k(node.pos)
    sp = ForwardDiff.derivative(totalsource,node.T)
    sc = totalsource(node.T) - node.T*sp
    an = []
    kn = 0.0
    del_n = 0.0
    for i = 1:1:length(node.neighbors)
        neighbornode = s.mesh.nodes[node.neighbors[i]]
        kN = s.k(neighbornode.pos)
        boundary = s.mesh.boundaries[node.boundaries[i]]
        pos_vector = neighbornode.pos - node.pos #Points from center of node to center of neighboring node
        del_n = norm(pos_vector)
        del_n_P = norm(boundary.pos - node.pos)#amount of del_n in the P region
        del_n_N = del_n - del_n_P
        kn = del_n / (del_n_P/kP + del_n_N/kN)
        pos_unit_vector = pos_vector / del_n
        V = s.velocity(node.pos)
        D = kn/del_n
        F = s.rho(node.pos)*dot(V,pos_unit_vector)
        P = F/D
        AP = 0
        if s.convectionscheme == "central"
            AP = 1 - 0.5*abs(P)
        elseif s.convectionscheme == "upwind"
            AP = 1
        elseif s.convectionscheme == "hybrid"
            AP = maximum([1 - 0.5*abs(P),0])
        elseif s.convectionscheme == "power"
            AP = maximum([(1 - 0.1*abs(P))^5,0])
        elseif s.convectionscheme == "exponential"
            AP = abs(P)/(exp(abs(P))-1)
        end
        #println(s.area(node.pos))
        push!(an,(D*AP+maximum([-F,0]))*s.area(node.pos))
        #NOTE: this assumes that  neighbor and boundary indecies correspond
        
    end
    #println(node.vol)
    ap = sum(an) - node.vol * sp
    b = node.vol *sc

    #overwrite the values if there is a BC
    if node.BC != 0
        BC = s.BCs[node.BC]
        if typeof(BC) <:flux
            ap = kn / del_n
            an = [kn / del_n]
            b = BC.value
        elseif typeof(BC) <:convective
            ap = kn / del_n + BC.value
            an = [kn / del_n]
            b = BC.value*BC.Tinf
        elseif typeof(BC) <:constanttemp
            ap = 1.0
            an=[0]
            b = BC.value
        else
            error("Invalid boundary condition type.")
        end
    end
    return ap,an,b
end


function set_up!(s::T) where T<:meshedoneDscene
    for i = 1:1:length(s.mesh.nodes)
        node = s.mesh.nodes[i]
        node.ap,node.an,node.b = find_a_b(s,i)
    end
end

function getA_b(s::T) where T<:meshedoneDscene
    n = length(s.mesh.nodes)
    A = zeros(n,n)
    b = zeros(n,1)
    for i = 1:1:n
        node =s.mesh.nodes[i]
        A[i,i] = -node.ap
        for j = 1:1:length(node.neighbors)
            neighborindex = node.neighbors[j]
            A[i,neighborindex] = node.an[j]
        end
        b[i,1] = - node.b
    end
    return A, b
end

function itterate_solve_steady!(s::T1;reltol =1e-4) where T1<:meshedoneDscene
    newTemps = fill(1E10,length(s.mesh.nodes)) #10^10 K is pretty unbelievable
    oldTemps = deepcopy(newTemps)
    error = 100
    itters = 0
    while error >= reltol
        set_up!(s)
        A,b = getA_b(s)
        global Amatrix = A
        global bee = b
        newTemps = tdma_solve(A,b)
        #newTemps = A \ b
        for i = 1:1:length(s.mesh.nodes)
            s.mesh.nodes[i].T = newTemps[i]
        end
        #NOTE: changed way to calculate error with HW3
        error = maximum(@. ((newTemps-oldTemps)^2)^0.5)#sum(@. ((newTemps - oldTemps)^2)/(oldTemps^2))^.5
        oldTemps = deepcopy(newTemps)
        itters = itters +1
        #println(itters) #Usually converges in 6 itters or less
    end
end

function itterate_solve_GS!(s::T1,newTemps;reltol =1e-6) where T1<:meshedoneDscene
    #newTemps = fill(1E10,length(s.mesh.nodes)) #10^10 K is pretty 
    for i = 1:1:length(s.mesh.nodes)
        s.mesh.nodes[i].T = newTemps[i]
    end
    oldTemps = deepcopy(newTemps)
    error = 100
    itters = 0
    while error >= reltol
        #set_up!(s)
        #A,b = getA_b(s)
        #newTemps = tdma_solve(A,b)
        #newTemps = A \ b
        for i = 1:1:length(s.mesh.nodes)
            neighbortemps = []
            for j = 1:1:length(s.mesh.nodes[i].neighbors)
                neighbornode = s.mesh.nodes[s.mesh.nodes[i].neighbors[j]]
                push!(neighbortemps,neighbornode.T)
            end
            ap,an,b = find_a_b(s,i)
            #println(neighbortemps)
            
            calc = (dot(an,neighbortemps) + b)/ap
            newTemps[i] = newTemps[i] + 1.0 * (calc - newTemps[i])
            s.mesh.nodes[i].T = newTemps[i]
        end
        #NOTE: changed way to calculate error with HW3
        #error = maximum(@. abs(newTemps-oldTemps))#1a
        dx = 0.003/(s.meshingsettings.ncells-2)
        k = s.k(0)
        qnew = 2*k*(newTemps[1]-newTemps[2])/dx
        qold = 2*k*(oldTemps[1]-oldTemps[2])/dx
        error = abs(qnew-qold)
        #println(error2)
        oldTemps = deepcopy(newTemps)
        itters = itters +1
        println(itters) #Usually converges in 6 itters or less
    end
end

function convergencestudy_steady(s::T1;reltol =1e-4) where T1<:oneDscene
    intsteps = 1E6
    newTemps = fill(1E10,3) #10^10 K is pretty unbelievable
    L = s.length
    newfit = linear_interpolation([0,L/2,L],newTemps)
    oldTemps = deepcopy(newTemps)
    errors = [100.0]
    itters = [0]
    CVsused = [0]
    nCV = 1
    while errors[end] >= reltol
        s.meshingsettings.ncells = nCV +2
        sm = mesh(s)
        set_up!(sm)
        A,b = getA_b(sm)
        newTemps = tdma_solve(A,b)#newTemps = A \ b
        newxs = []
        newTs = []
        for i = 1:1:length(sm.mesh.nodes)
            sm.mesh.nodes[i].T = newTemps[i]
            push!(newTs,newTemps[i])
            push!(newxs,sm.mesh.nodes[i].pos[1])
        end
        oldfit = deepcopy(newfit)
        newfit = linear_interpolation(newxs, newTs)
        function differencesquared(x)
            return (newfit(x) - oldfit(x))^2
        end
        Tavg = integrate(oldfit,0,L,intsteps)/(L)
        if nCV > 2
            intsteps =10*nCV*(nCV-1)
        else
            intsteps = 10
        end
        push!(errors,integrate(differencesquared,0,L,intsteps)^0.5/(Tavg*L))
        oldTemps = deepcopy(newTemps)
        push!(itters,itters[end] +1)
        push!(CVsused,length(sm.mesh.nodes)-2)
        nCV = nCV + 1
        #println(errors[end]) 
        println(newTemps[end])
    end
    
    return dedupe_and_correlate(CVsused[2:end],errors[2:end],xlabel = "Number of Control Volumes", ylabel = "Error",yaxis = :log)
end

function convergencestudy_knownsolution(s,f;reltol =1e-4)
    errors = [100.0]
    itters = [0]
    CVsused = [0]
    nCV = 5
    while errors[end] >= reltol
        s.meshingsettings.ncells = nCV +2
        sm = mesh(s)
        set_up!(sm)
        A,b = getA_b(sm)
        newTemps = tdma_solve(A,b)#newTemps = A \ b
        xs = []
        for i = 1:1:length(sm.mesh.nodes)
            push!(xs,sm.mesh.nodes[i].pos[1])
        end

        error = 0.0
        for i = 1:1:length(newTemps)
            if abs((newTemps[i])-f(xs[i])) > error
                error = abs((newTemps[i])-f(xs[i]))
            end
        end

        push!(errors,error)
        push!(itters,itters[end] +1)
        push!(CVsused,length(sm.mesh.nodes)-2)
        nCV = nCV + 1
        println(errors[end]) 
    end
    p = dedupe_and_correlate(CVsused[2:end],errors[2:end],xlabel = "Number of Control Volumes", ylabel = "Error",yaxis = :linear)
    return p
end

function timestep!(scene,dt)#TODO: only the fully implicit method is shown, and only for a linear 1D problem
    dx = norm(scene.mesh.nodes[2].pos-scene.mesh.nodes[1].pos)
    T = record_T(scene)
    Tstar = record_T(scene)
    Told = record_T(scene)
    n =length(scene.mesh.nodes)
    A = zeros(n,n)
    b = zeros(n)
    k = scene.k(0.0)
    rho = scene.rho(0.0)
    cp = scene.cp(0.0)
    itterflag = true
    count = 0
    Trecord = [T]
    while (maximum(@. abs(T- Trecord[end])) > 10E-6) || itterflag
        count = count +1
        #println(count)
        itterflag = false
        push!(Trecord,T)
        
        for i = 1:1:n
            totalsource = scene.sources[1].func
            sp = ForwardDiff.derivative(totalsource,T[i])
            sc = totalsource(T[i]) - T[i]*sp
            if i == 1
                A[1,1] = 1.0
                b[1] = scene.BCs[1].value
            elseif i < n
                A[i,i-1] = -k*dt/dx
                A[i,i] = rho*cp*dx -sp*dx*dt + 2*k*dt/dx
                A[i,i+1] = -k*dt/dx
                b[i] = sc*dx*dt + rho*cp*dx*Told[i]
            else
                A[i,i-1] = -2*k*dt/dx
                A[i,i] = rho*cp*dx -sp*dx*dt + 2*k*dt/dx
                b[i] = sc*dx*dt + rho*cp*dx*Told[i]
            end
        end
        global Ay = A
        global Be = b
        T = (A)\(b)
        #@assert T != Tstar
    end

    for i = 1:1:n
        scene.mesh.nodes[i].T = T[i]
    end
end

function unsteady_solve!(scene,timesteps)
    Thist = []
    for i = 1:1:(length(timesteps)-1)
        push!(Thist,record_T(scene))
        dt = timesteps[i+1]-timesteps[i]
        timestep!(scene,dt)
    end
    push!(Thist,record_T(scene))
    return Thist
end

function convergencestudy_grid(s,CVrange)
    CVnums = collect(CVrange[1]:1:CVrange[2])
    endvals = []
    for i = 1:1:length(CVnums)
        nCV = CVnums[i]
        s.meshingsettings.ncells = nCV +2
        sm = mesh(s)
        set_up!(sm)
        A,b = getA_b(sm)
        newTemps = tdma_solve(A,b)#newTemps = A \ b
        push!(endvals,newTemps[end])
    end
    p = dedupe_and_correlate(CVnums[1:end],endvals[1:end],xlabel = "Number of Control Volumes", ylabel = "Tip Temperature (K)",yaxis = :linear)
    return p
end

function convergence_study_time(s,trange,steprange)
    stepnums = collect(steprange[1]:steprange[1]:steprange[2])
    scene = mesh(s)
    endvals = []
    for i = 1:1:length(stepnums)
        timesteps = collect(trange[1]:(1/stepnums[i]):1)
        Thistory = unsteady_solve!(scene,timesteps)
        push!(endvals,Thistory[end][end])
        stepnum = stepnums[i]
        println("Evaluated for $stepnum timesteps.")
    end
    p = dedupe_and_correlate(stepnums[1:end],endvals[1:end],xlabel = "Number of Timesteps", ylabel = "Tip Temperature (K) at t=1s",yaxis = :linear)
    return p
end