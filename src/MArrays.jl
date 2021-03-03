""" MArrays.jl contains the functions and methods for the evolution of initial
amplitudes in a one-dimensional masses array with harmonic coupling. Here, we
study different coupling laws, their evolution on time, and their energy content """

""" Circular arrangement with equal strength coupling """

"Amplitude and argument functions"
Ac(N) = 1/(N + 1)
ωc(m, N) = (π * m) * Ac(N)

"Kernel function"
kec(t, k, m, n, l, N) = cos(2 * t * (√k) * sin(ωc(m, N))) * cos(2 * ωc(m, N) * (n - l))

"qc is the main function who returns the amplitude in the _n_ position at the
_t_ time"
function qc(k, n, t, N, V::Vector)
    u = zeros(Float64,(N,N))
    for m in 0:(N-1)
        for l in 0:(N-1)
            u[l + 1, m + 1] = kec(t, k, m + 1, n, l + 1, N) * V[l + 1]
        end
    end
    return Ac(N) * sum(u)
end
