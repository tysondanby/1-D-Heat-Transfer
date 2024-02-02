function linear(x)
    return x
end

function begincosinespacing(x)
    return 1 - cos(0.5*pi*x)
end

function endcosinespacing(x)
    return cos(0.5*pi*(x-1))
end

function bicosinespacing(x) #increased meshing on both ends, less in middle.
    return .5 - .5*cos(pi*x)
end