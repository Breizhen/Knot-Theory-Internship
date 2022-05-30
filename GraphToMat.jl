using Graphs, SimpleWeightedGraphs

g = SimpleWeightedGraph(5)

sources = [1, 1, 1, 2, 3, 3, 4, 4]
dest    = [2, 3, 3, 5, 4, 5, 5, 5]
weights = -ones(length(sources))
g = SimpleWeightedGraph(sources, dest, weights; combine = +)

function graph2mat(g::SimpleWeightedGraph)
    if g == SimpleWeightedGraph(1)
        return Array{Int64}(undef, 0, 0)
    else
        cpt::Int64 = 0
        m::Matrix{Int64} = Matrix{Int64}(undef, nv(g), nv(g))
        for i in 1:(nv(g)-1)                #On each line...
            for j in 1:(i-1)                #...for each element before the diagonal element
                cpt += m[i,j]               #......We sum up the elements in the counter
            end
            for j in (i+1):nv(g)            #...for each element after the diagonal element
                m[i,j] = -g.weights[i,j]
                m[j,i] = m[i,j]
                cpt += m[i,j]
            end
            m[i,i] = -cpt
            cpt = 0
        end
        m[nv(g), nv(g)] = -sum(m[nv(g), 1:(end-1)])
        return m[2:end, 2:end]
    end
end
