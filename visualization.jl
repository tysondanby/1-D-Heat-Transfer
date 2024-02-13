function basic_temp_plot(s)
    global Ts = []
    global xs = []
    for i = 1:1:length(s.mesh.nodes)
        push!(xs,s.mesh.nodes[i].pos)
        push!(Ts,s.mesh.nodes[i].T)
    end
    return plot(xs,Ts,xlabel = "Distance (m)", ylabel = "Temperature (K)")
end


function printtable(series,labels)
    if length(series) != length(labels)
        error("Number of series must match number of labels")
    end
    numberdecimals = 9
    series_strings = []
    max_string_length_in_series = zeros(length(series))
    for i = 1:1:length(series)
        push!(series_strings,vectorvaluestostring(series[i],numberdecimals))
        lengths = @. length(series_strings[i])
        push!(lengths,length(labels[i]))
        max_string_length_in_series[i] = maximum(lengths)
    end
    labeline = ""
    bars = ""
    for i =1:1:length(series)
        numchars = length(labels[i])
        room = max_string_length_in_series[i]
        spaces = room - numchars
        for k = 1:1:spaces
            labeline = labeline*" "
            bars = bars*"-"
        end
        labeline = labeline*(labels[i])
        for j = 1:1:length(labels[i])
            bars = bars*"-"
        end
        if i != length(series)
            labeline = labeline*"|"
            bars = bars*"+"
        end
    end
    lines = [labeline,bars]
    for i = 1:1:length(series[1])
        line = ""
        for j = 1:1:length(series)
            numchars = length(series_strings[j][i])
            room = max_string_length_in_series[j]
            spaces = room - numchars
            for k = 1:1:spaces
                line = line*" "
            end
            line = line*(series_strings[j][i])
            if j != length(series)
                line = line*"|"
            end

        end
        push!(lines,line)
    end
    for i =1:1:length(lines)
        println(lines[i])
    end
    #return lines
end