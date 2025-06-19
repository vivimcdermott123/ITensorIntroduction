using ITensors, ITensorMPS
using Plots

# Heisenberg Hamiltonian as MPO
function heisenberg_mpo(N; J=1.0, S="S=1/2")
    sites = siteinds(S, N)
    os = OpSum()
    for i in 1:N-1
        os += 0.5J, "S+", i, "S-", i+1
        os += 0.5J, "S-", i, "S+", i+1
        os += J, "Sz", i, "Sz", i+1
    end
    return MPO(os, sites), sites
end

# Run DMRG to find the ground state energy and wavefunction
function run_dmrg(H, sites; maxdim=100)
    init_state = [isodd(i) ? "Up" : "Dn" for i in 1:length(sites)]
    psi0 = MPS(sites, init_state)
    sweeps = Sweeps(10)
    setmaxdim!(sweeps, maxdim)
    energy, psi = dmrg(H, psi0, sweeps;)
    return energy, psi
end

# Calculate energy gap Δ(L) for given size L
function calculate_gap(L; J=1.0, maxdim=100)
    H, sites = heisenberg_mpo(L; J=J)
    E0, psi_gs = run_dmrg(H, sites; maxdim=maxdim)
    E1, psi_ex = dmrg(H, [psi_gs]; nsweeps=10, maxdim=maxdim, silent=true)
    return E1 - E0
end

# Correlation function ⟨Sᶻᵢ Sᶻⱼ⟩
function correlation(psi, sites, i, j)
    Sz_i = op(sites, "Sz", i)
    Sz_j = op(sites, "Sz", j)
    return inner(psi, Sz_i * Sz_j, psi)
end

# Entanglement entropy calculation
function entanglement_entropy(psi, site)
    orthogonalize!(psi, site)
    _, S, _ = svd(psi[site], (linkind(psi, site-1),))
    return -sum(s^2 * log(s^2) for s in diag(S))
end

# Magnetization profile ⟨Sz⟩
function magnetization_profile(psi, sites)
    [real(expect(psi, "Sz", i)) for i in 1:length(sites)]
end

# Plotting 1-5
function plot_all()
    ## 1. Energy Gap Δ(L)
    Ls = 4:2:16
    gaps = [calculate_gap(L) for L in Ls]
    plt1 = plot(Ls, gaps, marker=:o, xlabel="System Size L", ylabel="Energy Gap Δ",
        title="1. Energy Gap vs. System Size", legend=false)

    ## 2. Spin-Spin correlation
    L = 20
    H, sites = heisenberg_mpo(L)
    _, psi_gs = run_dmrg(H, sites)
    center = div(L, 2)
    corrs = [correlation(psi_gs, sites, center, j) for j in 1:L]
    plt2 = plot(1:L, corrs, xlabel="Site j", ylabel="⟨Sz_i Sz_j⟩",
        title="2. Spin-Spin Correlation", marker=:circle, legend=false)

    ## 3. Entanglement Entropy across bonds
    ent_entropy = [entanglement_entropy(psi_gs, i) for i in 2:L-1]
    plt3 = plot(2:L-1, ent_entropy, xlabel="Bond index",
        ylabel="Entanglement Entropy", title="3. Entanglement Entropy", marker=:diamond, legend=false)

    ## 4. Magnetization profile
    magnetization = magnetization_profile(psi_gs, sites)
    plt4 = bar(1:L, magnetization, xlabel="Site", ylabel="⟨Sz⟩",
        title="4. Magnetization Profile")

    ## 5. Excited vs Ground state local observables
    _, psi_ex = dmrg(H, [psi_gs]; nsweeps=10, maxdim=100, silent=true)
    magnetization_ex = magnetization_profile(psi_ex, sites)
    plt5 = plot(1:L, magnetization, label="Ground State", xlabel="Site",
        ylabel="⟨Sz⟩", title="5. Excited vs Ground State Magnetization", marker=:o)
    plot!(1:L, magnetization_ex, label="Excited State", marker=:square)

    # Layout all plots together
    plot(plt1, plt2, plt3, plt4, plt5, layout=(3,2), size=(1200,1000))
end

# Run plot_all to generate and display plots
plot_all()
