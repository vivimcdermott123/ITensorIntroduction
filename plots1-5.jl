using ITensors
using ITensorMPS
using Plots
# Construct the Heisenberg Hamiltonian as an MPO
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

# Run DMRG given an initial product state
function run_dmrg(H, psi0; maxdim=100, nsweeps=10)
    sweeps = Sweeps(nsweeps)
    setmaxdim!(sweeps, maxdim)
    energy, psi = dmrg(H, psi0, sweeps)
    return energy, psi
end

# Initialize state with fixed number of up spins (total Sz)
function initial_state_fixed_sz(sites, nup)
    N = length(sites)
    return [i <= nup ? "Up" : "Dn" for i in 1:N]
end

# Calculate energy gap Δ(L)
function calculate_gap_sz(L; J=1.0, maxdim=100)
    H, sites = heisenberg_mpo(L; J=J)
    psi0_gs = MPS(sites, initial_state_fixed_sz(sites, div(L,2)))
    psi0_ex = MPS(sites, initial_state_fixed_sz(sites, div(L,2)+1))
    E0, _ = run_dmrg(H, psi0_gs; maxdim=maxdim)
    E1, _ = run_dmrg(H, psi0_ex; maxdim=maxdim)
    return E1 - E0
end

# Proper calculation of full correlation matrix ⟨SzᵢSzⱼ⟩
function correlation_profile(psi, sites, center)
    corr_matrix = correlation_matrix(psi, "Sz", "Sz"; sites=sites, site_range=1:length(sites), ishermitian=true)
    return corr_matrix[center, :]
end

# Compute entanglement entropy across a bond
function entanglement_entropy(psi, b)
    orthogonalize!(psi, b)
    U, S, V = svd(psi[b], (linkind(psi, b-1),))
    s = diag(S)
    s = s[s .> 1e-12]
    return -sum(s.^2 .* log.(s.^2))
end

# Corrected magnetization profile calculation ⟨Sz⟩
function magnetization_profile(psi, sites)
    return [real(expect(psi, "Sz"; sites=sites, site_range=i)[1]) for i in 1:length(sites)]
end

# Main function to plot all analyses
function plot_all()
    # 1. Energy Gap Δ(L)
    Ls = 4:2:16
    gaps = [calculate_gap_sz(L) for L in Ls]
    plt1 = plot(Ls, gaps, marker=:o, xlabel="System Size L", ylabel="Energy Gap Δ",
        title="1. Energy Gap vs System Size", legend=false)

    # Setup single L for other analyses
    L = 10
    H, sites = heisenberg_mpo(L)
    psi0_gs = MPS(sites, initial_state_fixed_sz(sites, div(L,2)))
    _, psi_gs = run_dmrg(H, psi0_gs; maxdim=100)

    # 2. Spin-Spin Correlation ⟨Szᵢ Szⱼ⟩
    center = div(L,2)
    corrs = correlation_profile(psi_gs, sites, center)
    plt2 = plot(1:L, corrs, xlabel="Site j", ylabel="⟨Sz_center Sz_j⟩",
        title="2. Spin-Spin Correlation", marker=:circle, legend=false)

    # 3. Entanglement Entropy S(bond)
    ent_entropy = [entanglement_entropy(psi_gs, b) for b in 1:L-1]
    plt3 = plot(1:L-1, ent_entropy, xlabel="Bond Index", ylabel="Entanglement Entropy",
        title="3. Entanglement Entropy", marker=:diamond, legend=false)

    # 4. Magnetization Profile ⟨Sz⟩
    magnetization = magnetization_profile(psi_gs, sites)
    plt4 = bar(1:L, magnetization, xlabel="Site", ylabel="⟨Sz⟩",
        title="4. Magnetization Profile")

    # 5. Excited vs Ground State Magnetization
    psi0_ex = MPS(sites, initial_state_fixed_sz(sites, div(L,2)+1))
    _, psi_ex = run_dmrg(H, psi0_ex; maxdim=100)
    magnetization_ex = magnetization_profile(psi_ex, sites)
    plt5 = plot(1:L, magnetization, label="Ground State", xlabel="Site",
        ylabel="⟨Sz⟩", title="5. Magnetization: Ground vs Excited", marker=:o)
    plot!(plt5, 1:L, magnetization_ex, label="Excited State", marker=:square)

    # Arrange all plots
    plot(plt1, plt2, plt3, plt4, plt5, layout=(3,2), size=(1200,1000))
end

# Execute all plots
plot_all()
