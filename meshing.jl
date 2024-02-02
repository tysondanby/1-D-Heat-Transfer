
function mesh(s::T) where T <: oneDscene
    Tinit = s.Tsur
    boundaries = []
    nodes = []
    if s.meshingsettings.deploymentscheme == 'A'
        if length(s.layers) > 1
            @warn("Mesh deployment scheme A not recommended for composites. Please use scheme B.")
        end
        xs = collect(0:1/(s.meshingsettings.ncells-1):1)
        betweenxs = collect(.5/(s.meshingsettings.ncells-1):1/(s.meshingsettings.ncells-1):1)
        nodexs = @. s.meshingsettings.spacingfunc(xs) * s.length
        boundaryxs = @. s.meshingsettings.spacingfunc(betweenxs) * s.length
        for i = 1:1:length(boundaryxs)
            push!(boundaries,boundary(boundaryxs[i]))
        end
        push!(nodes,node(0,[2],[1],0.0,0.0,0.0,Tinit))
        for i = 2:1:(length(nodexs)-1)
            push!(nodes,node(nodexs[i],[i-1,i+1],[i-1,i],0.0,0.0,0.0,Tinit))
        end
        push!(nodes,node(nodexs[end],[length(nodexs)-1],[length(nodexs)-1],0.0,0.0,0.0,Tinit))


    else
        xs = collect(0:1/(s.meshingsettings.ncells-2):1)
        boundaryxs = @. s.meshingsettings.spacingfunc(xs) * s.length

        searchlength = 0 #INSERT ALL COMPOSITE BOUNDARIES
        layerindex = 1
        while (searchlength < s.length) && (layerindex <= length(s.layers))
            currentlayer = s.layers[layerindex]
            searchlength = searchlength + currentlayer.thickness
            if searchlength < s.length
                insert_and_dedup!(boundaryxs,searchlength)
                layerindex = layerindex + 1
            end
        end
        push!(boundaries,boundary(boundaryxs[1]))
        push!(nodes,node(0,[2],[1],0.0,0.0,0.0,Tinit))
        for i = 2:1:(length(boundaryxs)-1)
            push!(boundaries,boundary(boundaryxs[i]))
            push!(nodes,node((boundaryxs[i]+boundaryxs[i-1])/2,[i-1,i+1],[i-1,i],0.0,0.0,0.0,Tinit))
        end
        push!(boundaries,boundary(boundaryxs[end]))
        push!(nodes,node((boundaryxs[end]+boundaryxs[end-1])/2,[length(boundaryxs)-1,length(boundaryxs)+1],[length(boundaryxs)-1,length(boundaryxs)],0.0,0.0,0.0,Tinit))
        push!(nodes,node(boundaryxs[end],[length(boundaryxs)],[length(boundaryxs)],0.0,0.0,0.0,Tinit))
    end


    return meshedoneDscene(s.meshingsettings,s.length,s.area,s.sources,s.BCs,s.layers,s.Tsur,oneDmesh(nodes,boundaries),s.k)
end