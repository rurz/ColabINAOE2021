using ColabINAOE2021

export dqc, dql, dqk
export wijct, wijlt, wijkt

""" Energy exchanges """

# 1. ---------------------------------------------------------------------------

"""### Circular array energy exchange """

""" Redefine the equation qc() """
qc(k, n, t, N, V) = qc(k, n, t, N, V)

""" Define the time derivative of qc() """
function dqc(k, n, t, N, V::Vector)
    v = zeros(Float64, (N + 1, N + 1))
	for m in 0:N
		for l in 0:m
			v[l + 1, m + 1] = sin(2 * t * (√k) * sin(π * m / (N + 1))) * sin(π * m / (N + 1)) * cos(2 * π * m * (n - l) / (N + 1)) * V[l + 1]
		end
	end
	return ((-2 * (√k)) / (N + 1)) * sum(v)
end

""" Define the energy exchange rate. """
wijc(k, i, j, t, N, V) = qc(k, i, t, N, V) * dqc(k, j, t, N, V) - dqc(k, i, t, N, V) * qc(k, j, t, N, V)

""" Energy exchange vector between the elements _i_ and _j_ """
wijct(k, i, j, t, N, V) = [wijc(k, i, j, t[τ], N, V) for τ in 1:length(t)]

""" Calculate the energy exchange distribution between the element _i_ and all
others in the array """
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

function dql(k, n, t, N, V::Vector)
	v = zeros(Float64, (N + 1, N + 1))
	for i in 0:N
		for j in 0:N
			v[i + 1, j + 1] = (ColabINAOE2021.twocheb(j, n, i, N) / ColabINAOE2021.nc(j, N)) * sin(2 * t * (√k) * sin(ColabINAOE2021.ϕ(j + 1, N) / 2)) * sin(ColabINAOE2021.ϕ(j + 1, N) / 2) * V[i + 1]
		end
	end
	return (-2 * (√k)) * sum(v)
end

wijl(k, i, j, t, N, V) = ql(k, i, t, N, V) * dql(k, j, t, N, V) - dql(k, i, t, N, V) * ql(k, j, t, N, V)

wijlt(k, i, j, t, N, V) = [wijl(k, i, j, t[τ], N, V) for τ in 1:length(t)]

# 3. ---------------------------------------------------------------------------

"""### Kravchuk array energy exchange """

function dqk(n, t, N, V::Vector)
	v = zeros(Float64, (N + 1, N + 1))
	for m in 0:N
		for l in 0:N
			v[m + 1, l + 1] = (-√m) *ColabINAOE2021.k2(n, l, m, N) * sin(t * (√m)) * V[l + 1]
		end
	end
	return sum(v)
end

wijk(i, j, t, N, V) = qk(i, t, N, V) * dqk(j, t, N, V) - dqk(i, t, N, V) * qk(j, t, N, V)

wijkt(i, j, t, N, V) = [wijk(i, j, t[τ], N, V) for τ in 1:length(t)]
