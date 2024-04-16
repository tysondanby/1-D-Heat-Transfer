
function find_nodeindex_from_pos(s,nodes,pos)#finds which node is closest to the point "pos."
    #This can likely be sped up with a better search algorithm, but this one works for all cases, and is simple to implement.
    mindist = 1e10
    indexmindist = 0
    for i = 1:1:length(nodes)
        dist = norm(nodes[i].pos-pos)
        if dist < mindist
            indexmindist = i
            mindist = dist
        end
    end
    return indexmindist
end


function mesh(s::T) where T <: oneDscene
    #Tinit = s.Tsur
    boundaries = []
    nodes = []
    if s.meshingsettings.deploymentscheme == 'A'
        if length(s.layers) > 1
            @warn("Mesh deployment scheme A not recommended for composites. Please use scheme B.")
        end
        xs = collect(0:1/(s.meshingsettings.ncells-1):1)
        betweenxs = collect(.5/(s.meshingsettings.ncells-1):1/(s.meshingsettings.ncells-1):1)
        nodexs = fill([0.0,0.0,0.0],length(xs))
        for i = 1:1:length(xs)
            nodexs[i]=[s.meshingsettings.spacingfunc(xs[i]) * s.length,0,0]
        end
        boundaryposs = fill([0.0,0.0,0.0],length(betweenxs))
        for i = 1:1:length(betweenxs)
            boundaryposs[i]=[s.meshingsettings.spacingfunc(betweenxs[i]) * s.length,0,0]
            push!(boundaries,boundary(boundaryposs[i],s.area(boundaryposs[i])))
        end
        volume = s.area(nodexs[1]) * norm(boundaryposs[1]-nodexs[1])
        push!(nodes,node(nodexs[1],[2],[1],0.0,0.0,0.0,s.Tinit(nodexs[1]),volume,0))
        for i = 2:1:(length(nodexs)-1)
            volume = s.area(nodexs[i]) * norm(boundaryposs[i]-boundaryposs[i-1])
            push!(nodes,node(nodexs[i],[i-1,i+1],[i-1,i],0.0,0.0,0.0,s.Tinit(nodexs[i]),volume,0))
        end
        volume = s.area(nodexs[end]) * norm(boundaryposs[length(nodexs)-1]-nodexs[end])
        push!(nodes,node(nodexs[end],[length(nodexs)-1],[length(nodexs)-1],0.0,0.0,0.0,s.Tinit(nodexs[end]),volume,0))

    else
        xs = collect(0:1/(s.meshingsettings.ncells-2):1)
        boundaryposs = fill([0.0,0.0,0.0],length(xs))
        for i = 1:1:length(xs)
            boundaryposs[i] = [s.meshingsettings.spacingfunc(xs[i]) * s.length,0,0]
        end

        searchlength = 0 #INSERT ALL COMPOSITE BOUNDARIES
        layerindex = 1
        while (searchlength < s.length) && (layerindex <= length(s.layers))
            currentlayer = s.layers[layerindex]
            searchlength = searchlength + currentlayer.thickness
            if searchlength < s.length
                insert_and_dedup!(boundaryposs,searchlength)
                layerindex = layerindex + 1
            end
        end
        push!(boundaries,boundary(boundaryposs[1],s.area(boundaryposs[1])))
        push!(nodes,node([0,0,0],[2],[1],0.0,0.0,0.0,s.Tinit([0,0,0]),0.0,0))#Zero volume node
        for i = 2:1:(length(boundaryposs)-1)
            push!(boundaries,boundary(boundaryposs[i],s.area(boundaryposs[i])))
            volume = s.area((boundaryposs[i]+boundaryposs[i-1])/2) * norm(boundaryposs[i]-boundaryposs[i-1])
            push!(nodes,node((boundaryposs[i]+boundaryposs[i-1])/2,[i-1,i+1],[i-1,i],0.0,0.0,0.0,s.Tinit((boundaryposs[i]+boundaryposs[i-1])/2),volume,0))
        end
        push!(boundaries,boundary(boundaryposs[end],s.area(boundaryposs[end])))
        volume = s.area((boundaryposs[end]+boundaryposs[end-1])/2) * norm(boundaryposs[length(boundaryposs)]-boundaryposs[length(boundaryposs)-1])
        push!(nodes,node((boundaryposs[end]+boundaryposs[end-1])/2,[length(boundaryposs)-1,length(boundaryposs)+1],[length(boundaryposs)-1,length(boundaryposs)],0.0,0.0,0.0,s.Tinit((boundaryposs[end]+boundaryposs[end-1])/2),volume,0))
        push!(nodes,node(boundaryposs[end],[length(boundaryposs)],[length(boundaryposs)],0.0,0.0,0.0,s.Tinit(boundaryposs[end]),0.0,0))#Zero volume node
    end

    #Record the BCs
    for i = 1:1:length(s.BCs)
        BC = s.BCs[i]
        nodeindex = find_nodeindex_from_pos(s,nodes,BC.pos)
        nodes[nodeindex].BC = i
    end

    return meshedoneDscene(s.meshingsettings,s.convectionscheme,s.length,s.area,s.sources,s.BCs,s.layers,s.Tsur,oneDmesh(nodes,boundaries),s.k,s.velocity,s.rho,s.cp)
end