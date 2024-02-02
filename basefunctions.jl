insert_and_dedup!(v::Vector, x) = (splice!(v, searchsorted(v,x), [x]); v)

function integrate(f,lowb,upb,steps)#Midpoint method
    I = 0
    for i =1:2:(2*steps-1)
        x=lowb + i*(upb-lowb)/(steps*2)
        I = I + f(x)*(upb-lowb)/steps
    end
    return I
end