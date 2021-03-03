### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ b6b1b198-7bd8-11eb-10f5-a98dd0a59fcf
begin
	cd("..")
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 594227e6-7bd9-11eb-3ffb-ab8e6e7b4d95
begin
	using ColabINAOE2021
	using PyPlot
	using SparseArrays
	using LinearAlgebra
end

# ╔═╡ 2141105e-7bda-11eb-3519-1d5cff10d945
ion()

# ╔═╡ 6310aeb4-7bd9-11eb-0bd4-81a074081964
# Global Constants
begin
	k = 1
	N = 32
end;

# ╔═╡ 81e6e088-7bd9-11eb-17a6-03ef56dc026b
# Initial amplitude vector constructors

pint(q) = sparsevec(Dict(q + 1 => 1), N); # Single amplitude point at 0 ≤ q ≤ N

# ╔═╡ e3ceb6d0-7bda-11eb-115f-4584602f29df
# Initial amplitude vector

q0 = normalize(pint(16));

# ╔═╡ 203134d6-7bdb-11eb-2df1-374c190cc45b
begin 
	t_list = range(0,stop=60,length=200)
	lt = length(t_list)
end;

# ╔═╡ 77034402-7bdb-11eb-3b30-3dc93d022393
begin
	evolc = zeros(Float64,(N,lt))
	for n in 1:N
		for t in 1:lt
			evolc[n, t] = ColabINAOE2021.qc(k, n, t_list[t], N, Array(q0))
		end
	end
end

# ╔═╡ e30513a6-7bdb-11eb-0e31-cd6f16283698
begin
	fig_1 = figure(figsize = (10,5))
	ax_1 = gca()
	ax_1.imshow(evolc)
	ax_1.set_xlabel("Evolution time: t")
	ax_1.set_ylabel("Element index: n")
	ax_1.set_aspect("2")
	ax_1.set_xticks([0,50,100,150,200])
	ax_1.set_xticklabels([0,15,30,45,60])
	
	pcm = ax_1.get_children()[10]
	cb = colorbar(pcm,orientation="vertical",shrink=0.7,aspect=35,fraction=0.015)
	
	tight_layout()
	#savefig("notebooks/circ.png",dpi=150,transparent=true)
end

# ╔═╡ 07fccfdc-7bdc-11eb-1a00-a345a06f8517
sum(abs2.(evolc))

# ╔═╡ Cell order:
# ╠═b6b1b198-7bd8-11eb-10f5-a98dd0a59fcf
# ╠═594227e6-7bd9-11eb-3ffb-ab8e6e7b4d95
# ╠═2141105e-7bda-11eb-3519-1d5cff10d945
# ╠═6310aeb4-7bd9-11eb-0bd4-81a074081964
# ╠═81e6e088-7bd9-11eb-17a6-03ef56dc026b
# ╠═e3ceb6d0-7bda-11eb-115f-4584602f29df
# ╠═203134d6-7bdb-11eb-2df1-374c190cc45b
# ╠═77034402-7bdb-11eb-3b30-3dc93d022393
# ╠═e30513a6-7bdb-11eb-0e31-cd6f16283698
# ╠═07fccfdc-7bdc-11eb-1a00-a345a06f8517
