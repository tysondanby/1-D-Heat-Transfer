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

function set_up!(s::T) where T<:meshedoneDscene
    if s.meshingsettings.deploymentscheme == 'A'
        #Half node at begining
        @warn("Deployment Scheme A is not complete. use B.")
        #TODO: write a function for the half node's a,b
        for i = 2:1:(length(s.mesh.nodes)-1)
            node = s.mesh.nodes[i]
            #middle nodes
            node.ap,node.an,node.b = find_a_b_interior(s,i)
        end
        #half node at end

    else
        #zero-volume node at begining
        node = s.mesh.nodes[1]
        node.ap,node.an,node.b = find_a_b_zero_volume(s,1)
        for i = 2:1:(length(s.mesh.nodes)-1)
            node = s.mesh.nodes[i]
            #middle nodes
            node.ap,node.an,node.b = find_a_b_interior(s,i)
        end
        #zero-volume node at end
        node = s.mesh.nodes[end]
        node.ap,node.an,node.b = find_a_b_zero_volume(s,length(s.mesh.nodes))
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

function itterate_solve!(s::T1;reltol =1e-6) where T1<:meshedoneDscene
    newTemps = fill(1E10,length(s.mesh.nodes)) #10^10 K is pretty unbelievable
    oldTemps = deepcopy(newTemps)
    error = 100
    itters = 0
    while error >= reltol
        set_up!(s)
        A,b = getA_b(s)
        global Amatrix = A
        global bee = b
        newTemps = A \ b
        for i = 1:1:length(s.mesh.nodes)
            s.mesh.nodes[i].T = newTemps[i]
        end
        error = sum(@. ((newTemps - oldTemps)^2)/(oldTemps^2))^.5
        oldTemps = deepcopy(newTemps)
        itters = itters +1
        #println(itters) #Usually converges in 6 itters or less
    end
end

function convergencestudy(s::T1;reltol =1e-4) where T1<:oneDscene
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
        newTemps = A \ b
        newxs = []
        newTs = []
        for i = 1:1:length(sm.mesh.nodes)
            sm.mesh.nodes[i].T = newTemps[i]
            push!(newTs,newTemps[i])
            push!(newxs,sm.mesh.nodes[i].pos)
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
        println(errors[end]) 
    end
    p1 = plot(itters[3:end],errors[3:end], xlims = (0,itters[end]), xlabel = "Number of Iterations", ylabel = "Total Absolute Relative Error",yaxis=:log)
    p2 = plot(itters[3:end],CVsused[3:end], xlims = (0,itters[end]), ylims = (0,maximum(CVsused[3:end])), xlabel = "Number of Iterations",ylabel = "Number of Control Volumes")
    return p1,p2
end