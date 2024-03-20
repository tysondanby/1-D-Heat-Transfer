function basic_temp_plot(s)
    global Ts = []
    global xs = []
    for i = 1:1:length(s.mesh.nodes)
        push!(xs,s.mesh.nodes[i].pos[1])
        push!(Ts,s.mesh.nodes[i].T)
    end
    return plot(xs,Ts,xlabel = "Distance (m)", ylabel = "Temperature (K)")
end


function printtable(series,labels; decimals = 9,margin = 1)
    if length(series) != length(labels)
        error("Number of series must match number of labels")
    end
    numberdecimals = decimals
    series_strings = []
    leading_digits = []
    trailing_digits = []
    max_leading_digits = []
    max_trailing_digits = []
    title_biggerthan_numberby = []
    max_string_length_in_series = zeros(length(series))
    for i = 1:1:length(series)
        push!(series_strings,vectorvaluestostring(series[i],numberdecimals))
        lengths = @. length(series_strings[i])
        digits1 = @. findnumleadingchars(series_strings[i])
        digits2 = @. lengths - digits1 - 1
        push!(leading_digits, digits1)
        push!(trailing_digits, digits2)
        push!(max_leading_digits, maximum(digits1))
        push!(max_trailing_digits, maximum(digits2))
        push!(title_biggerthan_numberby,)
        push!(lengths,length(labels[i]))
        max_string_length_in_series[i] = maximum([lengths...,maximum(digits1)+1+maximum(digits2)])
        push!(title_biggerthan_numberby,maximum([max_string_length_in_series[i]-(maximum(digits1)+1+maximum(digits2)),0]))
    end
    labeline = ""
    bars = ""
    for i =1:1:length(series)
        numchars = length(labels[i])
        room = max_string_length_in_series[i]
        spaces = room - numchars
        for l = 1:1:margin
            labeline = labeline*" "
            bars = bars*"-"
        end
        labeline = labeline*(labels[i])
        for k = 1:1:spaces
            labeline = labeline*" "
            bars = bars*"-"
        end
        for j = 1:1:length(labels[i])
            bars = bars*"-"
        end
        for l = 1:1:margin
            labeline = labeline*" "
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
            #numchars = length(series_strings[j][i])
            #room = max_string_length_in_series[j]
            front_spaces = max_leading_digits[j] - leading_digits[j][i]
            rear_spaces = (max_trailing_digits[j] - trailing_digits[j][i])+ title_biggerthan_numberby[j]
            for l = 1:1:margin
                line = line*" "
            end
            for k = 1:1:front_spaces
                line = line*" "
            end
            line = line*(series_strings[j][i])
            for k = 1:1:rear_spaces
                line = line*" "
            end          
            for l = 1:1:margin
                line = line*" "
            end
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