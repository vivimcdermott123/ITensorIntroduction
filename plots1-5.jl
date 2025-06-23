julia
using ITensors, ITensorMPS, Plots
function plot_all()
           # 1. Energy Gap Δ(L)
           Ls = 4:2:16
           gaps = [calculate_gap_sz(L) for L in Ls]
           plt1 = plot(Ls, gaps, marker=:o, xlabel="System Size L", ylabel="Energy Gap Δ",
               title="1. Energy Gap vs System Size", legend=false)

           # Remaining plots use one chain length
           L = 100
           H, sites = heisenberg_mpo(L)
           nup_gs = div(L,2)
           psi0_gs = MPS(sites, initial_state_fixed_sz(sites, nup_gs))
           _, psi_gs = run_dmrg(H, psi0_gs; maxdim=100)

           # 2. Spin-Spin Correlation ⟨Sz_i Sz_j⟩
           center = div(L,2)
           corr_matrix = correlation_matrix(psi_gs, "Sz", "Sz"; sites=sites, site_range=1:L, ishermitian=true)
           corrs = corr_matrix[center, :]
           plt2 = plot(1:L, corrs, xlabel="Site j", ylabel="⟨Sz_center Sz_j⟩",
               title="2. Spin-Spin Correlation", marker=:circle, legend=false)

           # 3. Entanglement Entropy S(bond)
           ent_entropy = [entanglement_entropy(psi_gs, b) for b in 2:L-1]
           plt3 = plot(2:L-1, ent_entropy, xlabel="Bond Index", ylabel="Entanglement Entropy",
               title="3. Entanglement Entropy", marker=:diamond, legend=false)

           # 4. Magnetization Profile ⟨Sz⟩
           magnetization = magnetization_profile(psi_gs, sites)
           plt4 = bar(1:L, magnetization, xlabel="Site", ylabel="⟨Sz⟩",
               title="4. Magnetization Profile")

           # 5. Excited vs Ground State Magnetization Profile
           nup_ex = div(L,2) + 1
           psi0_ex = MPS(sites, initial_state_fixed_sz(sites, nup_ex))
           _, psi_ex = run_dmrg(H, psi0_ex; maxdim=100)
           magnetization_ex = magnetization_profile(psi_ex, sites)
           plt5 = plot(1:L, magnetization, label="Ground State", xlabel="Site",
               ylabel="⟨Sz⟩", title="5. Magnetization: Ground vs Excited", marker=:o)
           plot!(plt5, 1:L, magnetization_ex, label="Excited State", marker=:square)

           # Layout all plots in a grid
           plot(plt1, plt2, plt3, plt4, plt5, layout=(3,2), size=(1200,1000))
       end

       # Run the plotting function
plot_all (generic function with 1 method)
