#creatingnewfolder
mkdir ~/itensor_test
cd ~/itensor_test
#startjuliainfolder
julia
#activate/installITensor
import Pkg
Pkg.activate(".")           # activate local project
Pkg.add("ITensors")         # install ITensors cleanly
#tryloading
using ITensors
@which MPS

