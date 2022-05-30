#import Base.:+, *

function P(n::Int64) #returns a dictionary (a_i => b_i) s.t. it represents the polynomial P_n(A) = b_1*A^a_1 + ... + b_|n|*A^a_|n|
    m::Int64 = abs(n)
    sgn::Int = sign(n)
    polyn::Dict{Int64,Int64} = Dict{Int64,Int64}()
    if m!=0
        for i in 1:m
            polyn[sgn*(m - 4*i + 2)] = isodd(i) ? 1 : -1 #Ternary operator
        end
    end
    return polyn
end


##Transformations described in definition 1.2

#G'_ij
function trs1(i::Int64, j::Int64, G::Matrix{Int64})
    H = deepcopy(G)
    H[i,i] += H[i,j]
    H[j,j] += H[i,j]
    H[i,j] = 0
    H[j,i] = 0
    return H
end

#G''_ij
function trs2(i::Int64, j::Int64, H::Matrix{Int64})
    G = deepcopy(H)
    G[i,i] += G[j,j] + 2G[i,j]
    for k in (union(1:i-1, i+1:size(G,1))) #Works for i=1 and size(G,1)=1 => returns empty vector
        G[i,k] += G[j,k]
        G[k,i] += G[k,j]
    end
    return G[1:end .!= j, 1:end .!= j] #We delete the j-th row and column
end

#G'_i
function trs3(i::Int64, G::Matrix{Int64})
    return deepcopy(G[1:end .!= i, 1:end .!= i])
end

## μ function ##

function condition3(G::Matrix{Int64})
    i = 0
    test = false
    while (!test && i < size(G,1))
        i += 1
        H = G[i,:]
        test = test || (H[1:end .!= i] == zeros(length(H)-1))
    end
    if !test
        return -1
    else
        return i
    end
end

#Overloading sum for polynomials (Dict)
function Base.:+(A::Dict{Int64,Int64}, B::Dict{Int64,Int64})
    for i in keys(A)
        if haskey(B, i)
            B[i] += A[i]
        else
            B[i] = A[i]
        end
    end
    return B
end

#Overloading product for polynomials
function Base.:*(A::Dict{Int64, Int64}, B::Dict{Int64, Int64})
    C::Dict{Int64, Int64} = Dict{Int64, Int64}()
    for i in keys(A)
        for j in keys(B)
            if haskey(C, i+j)
                C[i+j] += A[i]*B[j]
            else
                C[i+j] = A[i]*B[j]
            end
        end
    end
    return C
end

#Removes the entries whose value is 0
function Base.:reduce(A::Dict{Int64, Int64})
    for i in keys(A)
        if A[i] == 0
            delete!(A,i)
        end
    end
    return A
end

function mu(G::Matrix{Int64})
    if isempty(G)
        return Dict(0 => 1)
    else
        i = condition3(G)
        if i != -1
            return reduce(((Dict(G[i,i] => 1))*(Dict(-2 => -1, 2 => -1)) + P(G[i,i]))*mu(trs3(i,G)))
        else
            vect = G[1,:]
            j = 2
            while vect[j] == 0
                j += 1
            end
            return reduce((Dict(-G[1,j] => 1))*mu(trs1(1,j,G)) + P(-G[1,j])*mu(trs2(1,j,G)))
        end
    end
end

#=
function mu(G::Matrix{Int64})
    return mu_rec(G, Dict{Int64,Int64}())
end
=#

#reduce(Dict(-2=>1)*mu(B) + Dict(4=>-1, 0=>1)*(Dict(-1=>1)*mu(F) + P(1)*mu(H)))
