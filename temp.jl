# Random assortment of functions for <kineticMC>
module kMCtemp

export boltzmann,ratecatalog,partialsums,deltat

# boltzmann constant, in eV K^{-1}
const kB = 8.6173324e-5 

function boltzmann(E::Float64,T::Float64)
    return exp(-E/(kB*T))
end

function createdislocation(N::Int64)
    # creates a straight dislocation line with <N> segments. For convenience,
    # we initialise it as an array of floating point values, though this may 
    # change later
    y = zeros(Int64,N)
end

function kinkwork(σ::Float64,a::Float64,b::Float64,h::Float64)
    return σ*a*b*h/2.0
end

function kinkevent(E::Float64,T::Float64,work::Float64,direction::Int64,v0::Float64)
    S = 3.0kB # entropy taken to be fixed
    free_energy = E - T*S - sign(direction)*work
    return v0*boltzmann(free_energy,T)
end

function ratecatalog(y::Array{Int64,1},h::Float64,
                     Jn::Array{Float64,1},Jm::Array{Float64,1})
    # initialise rate catalogue with length 2*<number of segments>
    # since r will frequently be zero, initialise to zero (check this, may
    # be better to assume that d1 and d2 are usually zero
    ####
    # NEED TO @assert Jn[i] > 0.0 and Jm[i] > 0.0
    N::Integer = length(y)                                 
    r = Array(Float64,2*N)
    
    # compute d1 = (y_{i-1} - y{i})/h, d2 = (y_{i+1} - y{i})
    d1 = Float64[float(y[i-1 > 0 ? i-1 : end] - y[i]) for i=1:N] 
    d2 = Float64[float(y[i+1 <= N ? i+1 : 1] - y[i]) for i=1:N]
    d = d1 + d2
    
    for i in 1:N
        if abs(d1[i]) < eps() && abs(d2[i]) < eps() # need to move function calls
            r[2i-1] = Jn[1]
            r[2i] = Jn[2]
        elseif d[i] > eps()
            r[2i-1] = Jm[1]*abs(d[i])
            r[2i] = 0.0 
        elseif d[i] < -eps()
            r[2i-1] = 0.0
            r[2i] = Jm[2]*abs(d[i])
        else # d1[i] = -d2[i] != 0.0
            r[2i] = 0.0
            r[2i-1] = 0.0
        end         
    end        
    return r
end

function partialsums(r::Array{Float64,1},R)
    ntransitions::Integer = length(r)
    sk = Array(Float64,ntransitions)
    
    # fill array with partial s[i] = Σ_{j=1}^{i}r_{j}/R
    for i=1:ntransitions
       sk[i] = sum(r[1:i])/R
    end
    return sk   
end

function totaltransitionrate(r::Array{Float64,1})
    R = sum(r)
end

function deltat(r::Array{Float64,1})
    # draws from the exponential distribution f(Dt) = R.exp(-R.Dt)
    R = totaltransitionrate(r)
    ξ1 = rand()
    Δt = -(1.0/R)*log(ξ1)
    return Δt
end

function chooseevent(r::Array{Float64,1})
    R = totaltransitionrate(r)
    sk = partialsums(r,R)
    eventk = 1
    
    # draw a random number from [0,1], and return k : s_{k-1} < ξ2 < s_{k}
    ξ2 = rand()
    if ξ2 < sk[1]
        nothing
    else
        for k=2:length(sk)
            if sk[k-1] < ξ2 <= sk[k]
                eventk = k
            end
        end
    end
    return eventk    
end

function executeevent(y::Array{Int64,1},k::Int64,h::Float64)
    p = mod(k,2)
    n::Int64 = (k+p)/2
    # note -> differs from Bulatov and Cai because Julia indexing starts at 1
    if isodd(k)
        # move segment n forward
        y[n] = y[n] + 1
    else # isodd(k) == true
        # move segment n backward
        y[n] = y[n] - 1
    end
    nothing # because y is mutable
end

function timestep(t::Float64,y::Array{Int64,1},rates::Array{Float64,1},
                                                              h::Float64)
    Δt = deltat(rates)
    k = chooseevent(rates)
    executeevent(y,k,h)
    return (t + Δt)
end

function runkineticMC(N::Integer,h::Float64,nsteps::Int64,Jn::Array{Float64,1},
                                                        Jm::Array{Float64,1})
    y = createdislocation(N)
    t = zeros(Float64,nsteps)
    meanloc = zeros(Float64,nsteps)
    for i=2:nsteps
        rates = ratecatalog(y,h,Jn,Jm)
        t[i] = timestep(t[i-1],y,rates,h)
        meanloc[i] = mean(y)
    end
    return t,y,meanloc
end

end # kMCtemp
