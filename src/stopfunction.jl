mutable struct StopFunction{T,S}
    oldvals::T
    counter::Int
    tol::S
    maxit::Int
end

function (st::StopFunction)(state)
    st.counter += 1
    st.counter > st.maxit && return true

    vals = state[2]
    diff = norm(vals - st.oldvals)
    diff <= st.tol && return true
    st.oldvals = vals

    return false
end