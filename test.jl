module kMCtest


function Whh()
    prefactor1 = μ/4π*(b1[1]*b2[1] + b1[2]*b1[2]/(1-υ))
    term1 = prefactor1*(f(x14,y13,rc)-f(x24,y13,rc)-f(x13,y13,rc)+f(x23,y13,rc))   
    term2 = μ/4π*(g(x14,y13,rc)-g(x24,y13,rc)-g(x13,y13,rc)+g(x23,y13,rc))
    return term1+term2 
end

function Whv()
    prefactor1 = μ/4π*(-b1[1]*b2[2] + (1-2υ)/(1-υ)*b1[2]*b2[1])
    term1 = prefactor1*(h(x13,y14,rc)-h(x23,y14,rc)-h(x13,y13,rc)+h(x23,y13,rc)
    prefactor2 = μ/4π*(b1[1]*b2[2]+2υ/(1-υ)*(b1[2]*b2[1]))
    term2 = prefactor2*(k(x13,y13,rc)-k(x23,y14,rc)-k(x13,y13,rc)+l(x23,y13,rc))
    return term1+term2
end

function f()
end

function g()
end

function h()
end

function k()
end

function E(W)
    # pass a system specific energy function W to E, to minimize clutter
end

end # kMCtest
