#Prepare initial state MPS
state = [iodd(n) ? "Up" : "Dn" for n=1:N]
psi0_i = MPS(sites,state)
#Do 10 sweeps of DMRG, gradually
#increasing the maximum MPS
# bond dimension
sweeps = Sweeps(10)
setmaxdim!(sweeps,10,20,100,200,400,800)
setcutoff!(sweeps,1E-8)

#Run the DMRG algorithm 
energy,psi0 = dmrg(H,psi0_i,sweeps)
