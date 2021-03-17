using SpecialPolynomials
using HypergeometricFunctions

export qc, ql, qk
export Qc, Ql, Qk

"MArrays.jl contains the functions and methods for the evolution of initial amplitudes in a one-dimensional masses array with harmonic couplings. Here, we study different coupling laws, their evolution on time, and their energy content (despicted at the EExchanges.jl script)"

# 1. ---------------------------------------------------------------------------

"### Circular arrangement of masses, _m_{j} = 1 ∀j_, for restitution constants _k_{j} = k, ∀j"

"`Ac(N)` gives the normalization term 1/(N + 1)"
Ac(N) = 1/(N + 1)
"`ωc(n, M)` gives an angular frequency who depends on position _m_"
ωc(m, N) = (π * m) * Ac(N)

"`kec` is the kernel function who provide the evolution in time. It is an oscillatory sinusoidal term who is called as `kec(t, k, m, n, l, N)`"
kec(t, k, m, n, l, N) = cos(2 * t * (√k) * sin(ωc(m, N))) * cos(2 * ωc(m, N) * (n - l))

"`qc` is the circular evolution function who integrates the kernel `kec` and their proper normalization terms `Ac` and `ωc`. For a constant _k_, it is called as `qc(k, n, t, N, V::Vector)`"
function qc(k, n, t, N, V::Vector)
    u = zeros(Float64,(N + 1, N + 1))
    for m in 0:N
        for l in 0:N
            u[l + 1, m + 1] = kec(t, k, m, n, l, N) * V[l + 1]
        end
    end
    return Ac(N) * sum(u)
end

"`Qc` is an ad-hoc function to obtain the matrix evolution from the function `qc`. It evaluates an initial _N + 1_ vector _V_ over a time span _t_. The resulting item is an _(N + 1) × `length(t)`_ matrix. It has middle performance. It is called `Qc(k, t, V::Vector)`"
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

"### Lineal arrays with specific conditions in the coupling matrix"

"##### Case m_{j} = 1, k_{j} = k, ∀j"

"We will use the generator method for obtain the second kind Chebyshev polynomial
given by the package SpecialPolynomials.jl"

"`gen(n) is the generator array for the Chebyshev's. It retuns a vector of _N + 1_ size with zeros everywhere, except at `gen(end) = 1`, who in turn serves as the generator for the n-th Chebyshev U polynomial."
function gen(n)
    gen = zeros(Int64, n + 1)
    gen[end] = 1
    return gen
end

"`chebU(x, n) gives the n-th Chebyshev polynomial U_{n}(x). It relies on the package SpecialPolynomials.jl"
function chebU(x, n)
    return SpecialPolynomials.ChebyshevU(gen(n))(x)
end

"`ϕ(j, N)` gives the j-th phase of the Chebyshev's U_{n}(x)"
ϕ(j, N) = (π * j) / (N + 2)
"`y(j, N)` gives the j-th zero of the Chebyshev's U_{n}(x)"
y(j, N) = cos(ϕ(j, N))

"`oscl` gives the j-th oscillatory evolution terms at time _t_, it is called as `oscl(k, t, j, N)`"
oscl(k, t, j, N) = cos(2 * t * (√k) * sin(ϕ(j + 1, N)/2))

"`nc(j, N)` gives a kind of _position dependent_ normalization, since it depends on the j-th term."
function nc(j, N)
    nc = zeros(N + 1)
    for l in 0:N
        nc[l + 1] = (chebU(y(j + 1, N), l))^2
    end
    return sum(nc)
end

"`twocheb` provides a compact form to package the multiplication of two mutually excludent Chebyshev's U_{n}(x)×U_{l}(x). It is called as `twocheb(j, l, n, N)`."
twocheb(j, l, n, N) = chebU(y(j + 1, N), l) * chebU(y(j + 1, N), n)

"`kel` is the composition of the `twocheb` and `oscl` functions. It gives the evolution kernel for the lineal array"
kel(t, k, n, l, j, N) = (twocheb(j, l, n, N)/nc(j, N)) * oscl(k, t, j, N)

"`ql` is the main function for the lineal array evolution, when evaluated for some initial _N + 1_ vector V, it gives a _(N + 1) × (N + 1)_ matrix who is a function time _t_ and position index _n_. It is called `ql(k, n, t, N, V::Vector)`"
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

"`Ql` is an ad-hoc function to obtain the matrix evolution from the function `ql`. It evaluates an initial _N + 1_ vector _V_ over a time span _t_. The resulting item is an _(N + 1) × `length(t)`_ matrix. It has poor performance. It is called `Ql(k, t, V::Vector)`"
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

"### Evolution in Kravchuk coupling. The coupling coefficients belongs to the basis of the discrete and finite harmonic oscillator of _su(2)_"

"`kp` gives the symmetric Kravchuk polynomial, it is evaluated as `kp(i, j, N)` for i,j ∈ [0,N]"
function kp(i, j, N)
    if i == j == N
        return HypergeometricFunctions._₂F₁general2(0, -i, -N, 2) # issue solved when i == j == N
    else
        return HypergeometricFunctions._₂F₁general2(-j, -i, -N, 2)
    end
end

"`ω(j, N)` gives a binomial normalization at the j-th term"
ω(j, N) = 2.0^(-N) * binomial(N, j)
"`h(i, N)` gives a binomial normalization at the i-th term"
h(i, N) = 1 / binomial(N, i)

"`oscc(t, n)` is the evolution kernel at time _t_. It is really simple"
oscc(t, n) = cos(t * sqrt(n))

"`kf` is the symmetric Kravchuk function, who is the composition of `ω` and `h`"
kf(i, j, N) = sqrt(ω(j, N) / h(i, N)) * kp(i, j, N)

"`kfm(N)` evaluates the Kravchuk matrix of `kf`. It returns a matrix of _(N + 1)×(N + 1)_ terms who columns are the functions `kf` at index _n_"
function kfm(N)
    km = zeros(Float64, (N + 1, N + 1))
    for i in 0:N
        for j in 0:N
            km[i + 1, j + 1] = kf(i, j, N)
        end
    end
    return km
end

"`k2` gives the two dimensional Kravchuk function, it is called `k2(m, l, n, N)`"
k2(m, l, n, N) = kf(m, n, N) * kf(l, n, N)

"`qk` is the main function for the Kravchuk array evolution, when evaluated for some initial _N + 1_ vector V, it gives a _(N + 1) × (N + 1)_ matrix who is a function time _t_ and position index _n_. It is called `qk(m, t, N, V::Vector)`"
function qk(m, t, N, V::Vector)
    u = zeros(Float64, (N + 1, N + 1))
    for l in  0:N
        for n in 0:N
            u[l + 1, n + 1] = k2(m, l, n, N) * oscc(t, n) * V[l + 1]
        end
    end
    return sum(u)
end

"`Qk` is an ad-hoc function to obtain the matrix evolution from the function `qk`. It evaluates an initial _N + 1_ vector _V_ over a time span _t_. The resulting item is an _(N + 1) × `length(t)`_ matrix. It has good performance. It is called `Qk(t, V::Vector)`"
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
