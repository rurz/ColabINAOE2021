using SpecialPolynomials
using HypergeometricFunctions

""" MArrays.jl contains the functions and methods for the evolution of initial
amplitudes in a one-dimensional masses array with harmonic coupling. Here, we
study different coupling laws, their evolution on time, and their energy content """

# 1. ---------------------------------------------------------------------------

"""### Circular arrangement with equal strength coupling: k_{j}\equiv j \forall j """

"* Amplitude and argument functions"
Ac(N) = 1/(N + 1)
ωc(m, N) = (π * m) * Ac(N)

"* Kernel function"
kec(t, k, m, n, l, N) = cos(2 * t * (√k) * sin(ωc(m, N))) * cos(2 * ωc(m, N) * (n - l))

"* qc is the main function who returns the amplitude in the _n_ position at the
_t_ time"
function qc(k, n, t, N, V::Vector)
    u = zeros(Float64,(N + 1, N + 1))
    for m in 0:N
        for l in 0:N
            u[l + 1, m + 1] = kec(t, k, m, n, l, N) * V[l + 1]
        end
    end
    return Ac(N) * sum(u)
end

"""* Qc is an ad-hoc function to evaluate the matrix evolution consisting on
a M-dimensional time array T and a N-dimensional initial amplitudes vector V_{0}.
So, we arrive to a matrix of N \times M. The function requires the value of the
restitution constant k, the time array T and the initial amplitudes vector V_{0} as
Qc(k, T:Array, V_{0}:Vector) """

function Qc(k, t, V::Vector)
    N = length(V)-1
    T = length(t)
    v = zeros(Float64,(N + 1,T))
    for l in 0:N
        for m in 1:T
            v[l + 1, m] = qc(k, l, t[m], N, V)
        end
    end
    return v
end

# The initial amplitudes vector V need to be of dimension N + 1 = dim,
# because the loops starts counting at 0.

# 2. ---------------------------------------------------------------------------

"""### Lineal array with specific conditions in the coupling matrix """

"""##### m_{j} = 1, k_{j} = k, \forall j """

""" We will use the generator method for obtain the second kind Chebyshev polynomial
given by the package SpecialPolynomials.jl """

"""* Generator array for the Chebyshev's """
function gen(n)
    gen = zeros(Int64, n + 1)
    gen[end] = 1
    return gen
end

"""* Chebyshev polynomial U_{n}(x) """
function chebU(x, n)
    return SpecialPolynomials.ChebyshevU(gen(n))(x)
end

""" Phase and zeros of the Chebyshev's """
ϕ(j, N) = (π * j) / (N + 2)
y(j, N) = cos(ϕ(j, N))

""" Oscillatory term """
oscl(k, t, j, N) = cos(2 * t * (√k) * sin(ϕ(j + 1, N)/2))

""" Position dependent normalization """
function nc(j, N)
    nc = zeros(N + 1)
    for l in 0:N
        nc[l + 1] = (chebU(y(j + 1, N), l))^2
    end
    return sum(nc)
end

""" Two-dimensional Chebyshev polymonial """
twocheb(j, l, n, N) = chebU(y(j + 1, N), l) * chebU(y(j + 1, N), n)

""" Kernel composition function """
kel(t, k, n, l, j, N) = (twocheb(j, l, n, N)/nc(j, N)) * oscl(k, t, j, N)

""" q_{l} is the main function for the lineal array evolution, when evaluated
    for some vector V of size N, it gives a N \times N matrix who is a function
    time t and position index n """
function ql(k, n, t, N, V::Vector)
    nr = [nc(j, N) for j in 0:N]
    tc(j, l, n) = twocheb(j, l, n, N)
    osc(k, t, j) = oscl(k, t, j, N)
    el = zeros(Float64,(N + 1, N + 1))
    for l in 0:N
        for j in 0:N
            el[l + 1, j + 1] = (tc(j, l, n) / nr[j + 1]) * osc(k, t, j) * V[l + 1]
        end
    end
    return sum(el)
end

""" Q_{l} is, again, an ad-hoc function for evalute the amplitude evolution
    function q_{l}(n, t) """
function Ql(k, t, V::Vector)
    dn = length(V)-1
    ell = zeros(Float64, (dn + 1, length(t)))
    for n in 0:dn
        for T in 1:length(t)
            ell[n + 1, T] = ql(k, n, t[T], dn, V)
        end
    end
    return ell
end

# 3. Kravchuk ------------------------------------------------------------------

""" Symmetric Kravchuk polynomial """
function kp(i, j, N)
    if i == j == N
        return HypergeometricFunctions._₂F₁general2(0, -i, -N, 2)
    else
        return HypergeometricFunctions._₂F₁general2(-j, -i, -N, 2)
    end
end

ω(j, N) = 2.0^(-N) * binomial(N, j)
h(i, N) = 1 / binomial(N, i)

oscc(t, n) = cos(t * sqrt(n))

""" Symmetric Kravchuk function """
kf(i, j, N) = sqrt(ω(j, N) / h(i, N)) * kp(i, j, N)

""" Kravchuk matrix """
function kfm(N)
    km = zeros(Float64, (N + 1, N + 1))
    for i in 0:N
        for j in 0:N
            km[i + 1, j + 1] = kf(i, j, N)
        end
    end
    return km
end

""" Two-dimensional Kravchuk function """
k2(m, l, n, N) = kf(m, n, N) * kf(l, n, N)

""" q_{k} function """
function qk(m, t, N, V::Vector)
    u = zeros(Float64, (N + 1, N + 1))
    for l in  0:N
        for n in 0:N
            u[l + 1, n + 1] = k2(m, l, n, N) * oscc(t, n) * V[l + 1]
        end
    end
    return sum(u)
end

""" Q_{k} is the ad-hoc function for straight-forward evaluation """
function Qk(t, V::Vector)
    dn = length(V)-1
    T = length(t)
    ek = zeros(Float64, (dn + 1, T))
    for m in 0:dn
        for τ in 1:T
            ek[m + 1, τ] = qk(m, t[τ], dn, V)
        end
    end
    return ek
end
