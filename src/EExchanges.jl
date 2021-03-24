using ColabINAOE2021

export dqc, dql, dqk
export wijct, wijlt, wijkt

"### `EExchanges.jl` gives the functions and method for the study of energy distribution and exchange within the evolution of initial amplitudes on masses arrays given in `MArrays.jl`"

# 1. ---------------------------------------------------------------------------

"### Circular array energy exchange "

"`dqc` is the derivative of the function `qc` for the evolution on circular masses arrays. It is called as `dqc(k, n, t, N, V::Vector)`"
function dqc(k, n, t, N, V::Vector)
    v = zeros(Float64, (N + 1, N + 1))
	for m in 0:N
		for l in 0:m
			v[l + 1, m + 1] = sin(2 * t * (√k) * sin(π * m / (N + 1))) * sin(π * m / (N + 1)) * cos(2 * π * m * (n - l) / (N + 1)) * V[l + 1]
		end
	end
	return ((-2 * (√k)) / (N + 1)) * sum(v)
end

"`wijc` is the energy exchange function between the position _i_ and _j_ at time _t_ of an initial amplitude vector _V_, such that, but not strictly, _i < j_. It is called `wijc(k, i, j, t, N, V::Vector)`."
wijc(k, i, j, t, N, V) = qc(k, i, t, N, V) * dqc(k, j, t, N, V) - dqc(k, i, t, N, V) * qc(k, j, t, N, V)

"`wijct` is an ad-hoc constructor (comprehension) to evalute the evolution of the energy exchange given a time span _t_ with an initial amplitudes vector _V_. It is called `wijct(k, i, j, t, N, V::Vector)`"
wijct(k, i, j, t, N, V) = [wijc(k, i, j, t[τ], N, V) for τ in 1:length(t)]

"`wijcm` is a function that calculates the energy exchange between the position _i_, and all other positions at time _t_. The item who return is a _N×`length(t)`_ matrix. It is called `wijcm(k, i, t, N, V::Vector)`"
function wijcm(k, i, t, N, V)
	v = zeros(Float64, (N, length(t)))
	for τ in 1:length(t)
			for j in 1:N
				v[j, τ] = wij(k, i, j, t[τ], N, V)
			end
		end
	return v
end

# 2. ---------------------------------------------------------------------------

"""### Lineal array energy exchange """

"`dql` is the derivative of the function `ql` for the evolution on lineal masses arrays. It is called as `dql(k, n, t, N, V::Vector)`"
function dql(k, n, t, N, V::Vector)
	v = zeros(Float64, (N + 1, N + 1))
	for i in 0:N
		for j in 0:N
			v[i + 1, j + 1] = (ColabINAOE2021.twocheb(j, n, i, N) / ColabINAOE2021.nc(j, N)) * sin(2 * t * (√k) * sin(ColabINAOE2021.ϕ(j + 1, N) / 2)) * sin(ColabINAOE2021.ϕ(j + 1, N) / 2) * V[i + 1]
		end
	end
	return (-2 * (√k)) * sum(v)
end

"`wijl` is the energy exchange function between the position _i_ and _j_ at time _t_ of an initial amplitude vector _V_, such that, but not strictly, _i < j_. It is called `wijl(k, i, j, t, N, V::Vector)`."
wijl(k, i, j, t, N, V) = ql(k, i, t, N, V) * dql(k, j, t, N, V) - dql(k, i, t, N, V) * ql(k, j, t, N, V)

"`wijlt` is an ad-hoc constructor (comprehension) to evalute the evolution of the energy exchange given a time span _t_ with an initial amplitudes vector _V_. It is called `wijlt(k, i, j, t, N, V::Vector)`"
wijlt(k, i, j, t, N, V) = [wijl(k, i, j, t[τ], N, V) for τ in 1:length(t)]

# 3. ---------------------------------------------------------------------------

"""### Kravchuk array energy exchange """

"`dqk` is the derivative of the function `ql` for the evolution on Kravchuk masses arrays. It is called as `dqk(n, t, N, V::Vector)`"
function dqk(n, t, N, V::Vector)
	v = zeros(Float64, (N + 1, N + 1))
	for m in 0:N
		for l in 0:N
			v[m + 1, l + 1] = (-√m) *ColabINAOE2021.k2(n, l, m, N) * sin(t * (√m)) * V[l + 1]
		end
	end
	return sum(v)
end

"`wijk` is the energy exchange function between the position _i_ and _j_ at time _t_ of an initial amplitude vector _V_, such that, but not strictly, _i < j_. It is called `wijk(i, j, t, N, V::Vector)`."
wijk(i, j, t, N, V) = qk(i, t, N, V) * dqk(j, t, N, V) - dqk(i, t, N, V) * qk(j, t, N, V)

"`wijkt` is an ad-hoc constructor (comprehension) to evalute the evolution of the energy exchange given a time span _t_ with an initial amplitudes vector _V_. It is called `wijkt(k, i, j, t, N, V::Vector)`"
wijkt(i, j, t, N, V) = [wijk(i, j, t[τ], N, V) for τ in 1:length(t)]
