### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ cb205c40-83a9-11eb-329c-e529ab101736
begin
	cd("..")
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 80030c2c-83ab-11eb-2a2e-39f1fcbcffc3
using ColabINAOE2021

# ╔═╡ 1c6fa178-83ad-11eb-154f-a11ce04426dd
using Plots

# ╔═╡ 5d4a918e-83b3-11eb-2445-f79669bef553
using LaTeXStrings

# ╔═╡ e90745a4-83ac-11eb-3c40-6f6d11de047e
using LinearAlgebra

# ╔═╡ 31a312be-83a8-11eb-365c-cd02c33f3b94
md" ## Energy exchange studies in the masses array evolution "

# ╔═╡ cdcb9f1a-83a8-11eb-14fe-67c783a8f2da
md"

In [1] we stipulate how to obtain the time evolution of initial amplitude arrays $\mathbf{Q}(0) = \sum_{l}q_{l}(0)\ket{l}$ under a sort of harmonic nearest-neighbour couplings given by coupling matrices $\mathbb{M}$. We can obtain closed-analytics solutions when the matrices $\mathbb{M}$ could be factorized as $\mathbb{A}\mathbb{B}\mathbb{A}^{\dagger}$ over a discrete and finite orthonormal space, like the Fock space of bras and kets.

-- Energy exchange: We define the energy exchange of one element _i_ to other element _j_ in an array as the function

$E_{i\rightarrow j} = k\left(q_{i}\dot{q}_{j} - \dot{q}_{i}q_{j}\right),$

for some restitution constant $k$. Because we have the explicit functions for three cases of arrays: one circular, and two lineals; we can obtain the energy exchange in such simple form, just obtaining the time-derivatives and composing the function.

To notice, is that the energy exchange depends on time, so we have specifically the energy exchange rate, that is, how much an element _j_ receive (or deliver) to an element _i_ in a timespan.

_____

[1] Urzúa, A. R., Ramos-Prieto, I., Soto-Eguibar, F., & Moya-Cessa, H. (2020). Dynamical analysis of mass–spring models using Lie algebraic methods. Physica A: Statistical Mechanics and Its Applications, 540, 123193. [doi](https://doi.org/10.1016/j.physa.2019.123193)
"

# ╔═╡ c00b788e-83ad-11eb-3456-1937a68cc057
pyplot()

# ╔═╡ 46a432d4-83ac-11eb-2410-93ca49a73b74
# Global Constants

begin
	N = 32
	dim = N + 1
	k = 1
end;

# ╔═╡ 7978a3de-83ac-11eb-3ec2-ef162c5ca306
t_list = range(0, stop = 200, length = 200)

# ╔═╡ b22933a6-83ac-11eb-36e0-53bda1f79066
begin
	v0 = zeros(dim)
	v0[1] = 1
	#v0[] = 1
	
	q0 = normalize(v0)
end;

# ╔═╡ 84833ef8-83b4-11eb-3b7d-0388d5b1a38e
md"### 1. Circular array "

# ╔═╡ f9f33a28-83b2-11eb-3b1c-f7db876a1e0a
plot([qc(1, 1, t, N, q0) for t in t_list], linewidth = 2, label = L"q_{1}(t)", xlabel = L"Time\; [t]", ylabel = L"Amplitude")

# ╔═╡ 0a5165a6-83af-11eb-3d16-91b1821ac83c
ec01 = wijct(1, 0, 1, t_list, N, q0);

# ╔═╡ 2e85e55c-83af-11eb-2b62-f7bcf8192c37
plot(ec01, linewidth = 2, label = L"E_{0\rightarrow 1}(t)", xlabel = L"Time\; [t]", ylabel = L"Energy")

# ╔═╡ 9616dc6a-83b4-11eb-3f55-3b3b13ae6a94
md" ### 2. Lineal array "

# ╔═╡ a103ac0a-83b4-11eb-0766-cdea31fc445e
plot([ql(1, 1, t, N, q0) for t in t_list], linewidth = 2, label = L"q_{1}(t)", xlabel = L"Time\; [t]", ylabel = L"Amplitude")

# ╔═╡ c646c242-83b4-11eb-009d-038d9bc9df7c
el01 = wijlt(1, 0, 1, t_list, N, q0);

# ╔═╡ ce72d0a0-83b4-11eb-11be-e38b4881c08d
plot(el01, linewidth = 2, label = L"E_{0\rightarrow 1}(t)", xlabel = L"Time\; [t]", ylabel = L"Energy")

# ╔═╡ 190f4c7e-83b5-11eb-2f75-7559d8264ede
md" ### 3. Kravchuk array "

# ╔═╡ 284ee190-83b5-11eb-0940-917c1c56f682
plot([qk(1, t, N, q0) for t in t_list], linewidth = 2, label = L"q_{1}(t)", xlabel = L"Time\; [t]", ylabel = L"Amplitude")

# ╔═╡ 3407b2f0-83b5-11eb-3003-6787a4694722
ek01 = wijkt(0, 1, t_list, N, q0);

# ╔═╡ a37ee8c8-83b6-11eb-32e4-51d7e3407b42
plot(ek01, linewidth = 2, label = L"E_{0\rightarrow 1}(t)", xlabel = L"Time\; [t]", ylabel = L"Energy")

# ╔═╡ e6099ea2-83b6-11eb-102c-45960096ab4e
[sum(ec01), sum(el01), sum(ek01)]

# ╔═╡ Cell order:
# ╟─31a312be-83a8-11eb-365c-cd02c33f3b94
# ╟─cdcb9f1a-83a8-11eb-14fe-67c783a8f2da
# ╠═cb205c40-83a9-11eb-329c-e529ab101736
# ╠═80030c2c-83ab-11eb-2a2e-39f1fcbcffc3
# ╠═1c6fa178-83ad-11eb-154f-a11ce04426dd
# ╠═5d4a918e-83b3-11eb-2445-f79669bef553
# ╠═c00b788e-83ad-11eb-3456-1937a68cc057
# ╠═e90745a4-83ac-11eb-3c40-6f6d11de047e
# ╠═46a432d4-83ac-11eb-2410-93ca49a73b74
# ╠═7978a3de-83ac-11eb-3ec2-ef162c5ca306
# ╠═b22933a6-83ac-11eb-36e0-53bda1f79066
# ╟─84833ef8-83b4-11eb-3b7d-0388d5b1a38e
# ╠═f9f33a28-83b2-11eb-3b1c-f7db876a1e0a
# ╠═0a5165a6-83af-11eb-3d16-91b1821ac83c
# ╠═2e85e55c-83af-11eb-2b62-f7bcf8192c37
# ╟─9616dc6a-83b4-11eb-3f55-3b3b13ae6a94
# ╠═a103ac0a-83b4-11eb-0766-cdea31fc445e
# ╠═c646c242-83b4-11eb-009d-038d9bc9df7c
# ╠═ce72d0a0-83b4-11eb-11be-e38b4881c08d
# ╟─190f4c7e-83b5-11eb-2f75-7559d8264ede
# ╠═284ee190-83b5-11eb-0940-917c1c56f682
# ╠═3407b2f0-83b5-11eb-3003-6787a4694722
# ╠═a37ee8c8-83b6-11eb-32e4-51d7e3407b42
# ╠═e6099ea2-83b6-11eb-102c-45960096ab4e
