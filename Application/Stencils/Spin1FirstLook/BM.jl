import SpiderWebModel as SW
using ThreadPinning
pinthreads(:cores)

S = SW.stencilConfig(zeros(32,32),1,boundaryCondition=:periodic)
ψG = SW.SimpleJastrowFunction(S)
CT = SW.ContinuousTimeMethod(0.1,w_avg_estimate = 0.2*length(S),Hxx = SW.Hxx_RK(0.1))
a = SW.startManyWalkerGFMC(S,CT,10*20,200,ψG);

@time a = SW.startManyWalkerGFMC(S,CT,100*20,200,ψG);
@time a = SW.startManyWalkerGFMC(S,CT,500*20,200,ψG);

nothing