insert_and_dedup!(v::Vector, x) = (splice!(v, searchsorted(v,x), [x]); v)

function integrate(f,lowb,upb,steps)#Midpoint method
    I = 0
    for i =1:2:(2*steps-1)
        x=lowb + i*(upb-lowb)/(steps*2)
        I = I + f(x)*(upb-lowb)/steps
    end
    return I
end

function richardson_extrapolate_results(outputs,nsCVs)
    #TODO: only works for spacing multiples of two
    order = round(log((outputs[2]-outputs[1])/(outputs[3]-outputs[2]))/log(2))
    est = outputs[3] + (outputs[2]-outputs[3])/(1-2^order)
    return order, est
end

function tdma_solve(A,be)
    d = 1*be
    a = [A[1,1]]
    b = [-A[1,2]]
    c = [0.0]
    P = [b[1]/a[1]]
    Q = [d[1]/a[1]]
    for i = 2:1:(length(A[:,1])-1)
        push!(a,A[i,i])
        push!(b,-A[i,i+1])
        push!(c,-A[i,i-1])
        push!(P,b[i]/(a[i]-(c[i]*P[i-1])))
        push!(Q,(d[i]+(c[i]*Q[i-1]))/(a[i]-(c[i]*P[i-1])))
    end
    push!(a,A[length(A[:,1]),length(A[:,1])])
    push!(b,0.0)
    push!(c,-A[length(A[:,1]),length(A[:,1])-1])
    push!(P,b[end]/(a[end]-(c[end]*P[end])))
    push!(Q,(d[end]+(c[end]*Q[end]))/(a[end]-(c[end]*P[end-1])))

    phi = zeros(length(d))
    phi[end] = Q[end]
    for i = 1:1:(length(d)-1)
        phi[end-i] = P[end-i]*phi[end-i+1]+Q[end-i]
    end
    return phi
end

function dedupe_and_correlate(X,Y;xlabel="",ylabel="",yaxis = yaxis)
    newX = []
    newY = []
    prevx = NaN
    for i = 1:1:length(X)
        if X[i] != prevx
            push!(newX,X[i])
            push!(newY,Y[i])
        end
        prevx = X[i]
    end
    return plot(newX,newY,xlabel=xlabel,ylabel=ylabel,yaxis = yaxis)
end

function vectorvaluestostring(v,nd)
    strings = []
    for i = 1:1:length(v)
        if typeof(v[i]) == String
            push!(strings,v[i])
        else
            temp = round(v[i]*(10^nd))/(10^nd)
            push!(strings,"$temp")
        end
    end
    return strings
end

function findnumleadingchars(str)
    index = 1
    for i = 1:1:length(str)
        if str[i] == '.'
            index = i
        end
    end
    return index - 1
end