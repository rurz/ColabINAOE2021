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

pint(q) = sparsevec(Dict(q + 1 => 1), N + 1); # Single amplitude point at 0 ≤ q ≤ N

# ╔═╡ e3ceb6d0-7bda-11eb-115f-4584602f29df
# Initial amplitude vector

q0 = normalize(pint(8) + pint(24));

# ╔═╡ 203134d6-7bdb-11eb-2df1-374c190cc45b
t_list = range(0, stop = 61, length = 200);

# ╔═╡ b2181d84-7fcf-11eb-2400-23cdf5b8fff7
md"### 1. Circular array"

# ╔═╡ b49a93f8-806b-11eb-2dc1-3f8e020d0f8c
evolc = ColabINAOE2021.Qc(1, t_list, Array(q0));

# ╔═╡ e30513a6-7bdb-11eb-0e31-cd6f16283698
begin
	fig_1 = figure(figsize = (10,5))
	ax_1 = gca()
	ax_1.imshow(evolc)
	ax_1.set_xlabel("Evolution time: t")
	ax_1.set_ylabel("Element index: n")
	ax_1.set_aspect("2")
	ax_1.set_xticks([0,33,66,99,132,165,200])
	ax_1.set_xticklabels([0,10,20,30,40,50,60])

	pcm1 = ax_1.get_children()[10]
	cb1 = colorbar(pcm1,orientation="vertical",shrink=0.7,aspect=35,fraction=0.015)

	tight_layout()
	#savefig("notebooks/circ.png",dpi=150,transparent=true)
end

# ╔═╡ 07fccfdc-7bdc-11eb-1a00-a345a06f8517
sum(abs2.(evolc))

# ╔═╡ ac233d02-7fcf-11eb-1640-c172d67a0cfc
md"### 2. Lineal array"

# ╔═╡ cbc567c2-8095-11eb-0996-05753de896af
evoll = ColabINAOE2021.Ql(1, t_list, Array(q0));

# ╔═╡ fb886f34-8095-11eb-0b4f-e58f8dcb8964
begin
	fig_2 = figure(figsize = (10,5))
	ax_2 = gca()
	ax_2.imshow(evoll)
	ax_2.set_xlabel("Evolution time: t")
	ax_2.set_ylabel("Element index: n")
	ax_2.set_aspect("2")
	ax_2.set_xticks([0,33,66,99,132,165,200])
	ax_2.set_xticklabels([0,10,20,30,40,50,60])

	pcm2 = ax_2.get_children()[10]
	cb2 = colorbar(pcm2,orientation="vertical",shrink=0.7,aspect=35,fraction=0.015)

	tight_layout()
	#savefig("notebooks/lin.png",dpi=150,transparent=true)
end

# ╔═╡ 64040048-8096-11eb-2c6f-ed400b9dbf42
md"### 3. Kravchuk array "

# ╔═╡ 6ec1d92e-8096-11eb-3980-7107878eecb7
evolk = ColabINAOE2021.Qk(t_list, Array(q0));

# ╔═╡ 7f1a6296-8096-11eb-12b1-a14e3a7eb885
begin
	fig_3 = figure(figsize = (10,5))
	ax_3 = gca()
	ax_3.imshow(evolk)
	ax_3.set_xlabel("Evolution time: t")
	ax_3.set_ylabel("Element index: n")
	ax_3.set_aspect("2")
	ax_3.set_xticks([0,33,66,99,132,165,200])
	ax_3.set_xticklabels([0,10,20,30,40,50,60])

	pcm3 = ax_2.get_children()[10]
	cb3 = colorbar(pcm3,orientation="vertical",shrink=0.7,aspect=35,fraction=0.015)

	tight_layout()
	#savefig("notebooks/krav.png",dpi=150,transparent=true)
end

# ╔═╡ Cell order:
# ╠═b6b1b198-7bd8-11eb-10f5-a98dd0a59fcf
# ╠═594227e6-7bd9-11eb-3ffb-ab8e6e7b4d95
# ╠═2141105e-7bda-11eb-3519-1d5cff10d945
# ╠═6310aeb4-7bd9-11eb-0bd4-81a074081964
# ╠═81e6e088-7bd9-11eb-17a6-03ef56dc026b
# ╠═e3ceb6d0-7bda-11eb-115f-4584602f29df
# ╠═203134d6-7bdb-11eb-2df1-374c190cc45b
# ╟─b2181d84-7fcf-11eb-2400-23cdf5b8fff7
# ╠═b49a93f8-806b-11eb-2dc1-3f8e020d0f8c
# ╠═e30513a6-7bdb-11eb-0e31-cd6f16283698
# ╠═07fccfdc-7bdc-11eb-1a00-a345a06f8517
# ╟─ac233d02-7fcf-11eb-1640-c172d67a0cfc
# ╠═cbc567c2-8095-11eb-0996-05753de896af
# ╠═fb886f34-8095-11eb-0b4f-e58f8dcb8964
# ╟─64040048-8096-11eb-2c6f-ed400b9dbf42
# ╠═6ec1d92e-8096-11eb-3980-7107878eecb7
# ╠═7f1a6296-8096-11eb-12b1-a14e3a7eb885
