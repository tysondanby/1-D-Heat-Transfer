abstract type Scene end
abstract type UnmeshedScene<:Scene end
abstract type MeshedScene<:Scene end

mutable struct oneDscene<:UnmeshedScene
    meshingsettings
    length
    area # Function cross-sectional area = a(x)
    sources
    BCs
    layers
    Tsur
    k
end

mutable struct meshedoneDscene<:MeshedScene
    meshingsettings
    length
    area # Function cross-sectional area = a(x)
    sources
    BCs
    layers
    Tsur
    mesh
    k #function of position giving thermal conductivity
end

struct meshingsettings
    deploymentscheme#A or B
    ncells#Number of cells
    spacingfunc# one to one function such that f(0) = 0 and f(1) = 1
end

struct source
    range
    func #function of node temperature
end

mutable struct node
    pos
    neighbors #vector of node index numbers
    boundaries #vector of edge index numbers
    ap
    an #vector. corresponds to entries in .neighbors
    b
    T
end

struct boundary #this will need to be quite different for 3D
    pos
end

abstract type BC end

struct flux <: BC
    value#flux per unit area
end

struct convective <: BC
    value#h
    Tinf
end

struct constanttemp <: BC
    value
end

mutable struct oneDmesh
    nodes
    boundaries
end

struct layer
    materialname
    thickness
end

function dummy()
    return
end