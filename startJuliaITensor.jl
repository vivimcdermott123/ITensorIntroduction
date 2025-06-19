rm -rf ~/.julia/environments
rm -rf ~/.julia/packages/ITensors

#relaunchingjulia 
julia
import Pkg
Pkg.add("ITensors")
#checkinstallation
using ITensors
@which MPS
