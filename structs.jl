abstract type Scene end
abstract type UnmeshedScene<:Scene end
abstract type MeshedScene<:Scene end

mutable struct oneDscene<:UnmeshedScene
    meshingsettings
    convectionscheme
    length
    area # Function cross-sectional area = a(x)
    sources
    BCs
    layers
    Tsur
    k
    velocity
    rho
    cp
    Tinit
end

mutable struct meshedoneDscene<:MeshedScene
    meshingsettings
    convectionscheme
    length
    area # Function cross-sectional area = a(x)
    sources
    BCs
    layers
    Tsur
    mesh
    k #function of position giving thermal conductivity
    velocity #function of position giving velocity
    rho #function of position giving density
    cp
end

mutable struct meshingsettings
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
    #NOTE: this assumes that  neighbor and boundary indecies correspond
    ap
    an #vector. corresponds to entries in .neighbors
    b
    T
    vol #volume
    BC #associated boundary condition (index). This is zero if no boundary condition is associated.
end

struct boundary #this will need to be quite different for 3D
    pos
    area 
end

abstract type BC end

struct flux <: BC
    pos
    value#flux per unit area
end

struct convective <: BC
    pos
    value#h
    Tinf
end

struct constanttemp <: BC
    pos
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