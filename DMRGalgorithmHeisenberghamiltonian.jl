julia
using ITensors, ITensorMPS

# Define a function that builds the Heisenberg spin-1/2 Hamiltonian as an MPO
function heisenberg_mpo(N)
  sites = siteinds("S=1/2", N)  # Create a vector of N spin-1/2 sites (physical indices)
  os = OpSum()                  # Initialize a symbolic operator sum

  # Add terms for nearest-neighbor Heisenberg interactions
  for i = 1:N-1
    os += "Sz", i, "Sz", i+1       # Sz_i * Sz_{i+1}
    os += 1/2, "S+", i, "S-", i+1  # (1/2) * S+_i * S-_i+1
    os += 1/2, "S-", i, "S+", i+1  # (1/2) * S-_i * S+_i+1
  end

  H = MPO(os, sites)  # Convert symbolic operator sum to an MPO using the site indices
  return H, sites     # Return both the MPO and the sites so we can use them later
end

H, sites = heisenberg_mpo(100)  # Build the Heisenberg MPO for a 100-site chain

# ---------------- Initial State ----------------

N = length(sites)  # Number of sites in the system (should be 100)

# Build a "Neel state" alternating up and down spins: ↑↓↑↓...
state = [isodd(n) ? "Up" : "Dn" for n in 1:N]

# Construct an MPS from the product state
psi0_i = MPS(sites, state)

# ---------------- DMRG Parameters ----------------

sweeps = Sweeps(10)  # Create a Sweeps object for 10 DMRG sweeps

# Set maximum MPS bond dimensions for each sweep
setmaxdim!(sweeps, 10, 20, 100, 200, 400, 800)

# Set truncation cutoff (controls accuracy vs. speed)
setcutoff!(sweeps, 1e-8)

# ---------------- Run DMRG ----------------

# Run the DMRG algorithm to minimize ⟨psi|H|psi⟩
# Starts from psi0_i, returns energy and optimized MPS
energy, psi0 = dmrg(H, psi0_i, sweeps)

# You can now inspect the ground state energy:
println("Ground state energy: ", energy)

