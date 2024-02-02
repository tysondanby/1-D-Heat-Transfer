function basic_temp_plot(s)
    global Ts = []
    global xs = []
    for i = 1:1:length(s.mesh.nodes)
        push!(xs,s.mesh.nodes[i].pos)
        push!(Ts,s.mesh.nodes[i].T)
    end
    return plot(xs,Ts,xlabel = "Distance (m)", ylabel = "Temperature (K)")
end