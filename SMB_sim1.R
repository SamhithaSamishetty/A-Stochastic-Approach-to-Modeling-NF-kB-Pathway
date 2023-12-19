library(ggplot2)
sim_name <- "SMB Final Project"
final_time <- 1000
# rate constants
params <- c(ka=1, ki=0.1, k1=1, k1i=0.1, k4=1, k5=0.1, k2=0.1, k2i=1, k3=0.1, k3i=0.1)
# initial conditions of the system
initial_state <- c(IKKa=0, IKKn=1, NF_kBa = 0, NF_kBi=1, IkBaa = 0, A20a = 0, IkBai=1, A20i=1, RNAs=0)
# defining reactions
reactions <- list(
  reaction("ka*IKKn*(1/(1+A20a)", c(IKKn = -1, IKKa=+1)), # IKK activation
  reaction("ki*IKKa*A20a", c(IKKa = -1, IKKn=+1)), # IKK inactivation
  reaction("k1*NF_kBi*IKKa*(1/(1+IkBaa))", c(NF_kBi=-1, NF_kBa=+1)), # NFkB activation
  reaction("k1i*NF_kBa*IkBaa", c(NF_kBa=-1, NF_kBi=+1)), # NFkB inactivation
  reaction("k2*NF_kBa*IkBai", c(IkBai=-1, IkBaa=+1)), # IkBa activation
  reaction("k2i*IKKa*A20a*IkBaa", c(IkBaa=-1, IkBai=+1)), # IkBa inactivation
  reaction("k3*A20i*NF_kBa", c(A20i=-1,A20a=+1)), # A20 activation
  reaction("k3i*A20a*NF_kBi", c(A20a=-1, A20i=+1)), # A20 inactivation
  reaction("k4*NF_kBa", c(RNAs=+1)),
  reaction("k5*RNAs", c(RNAs=-1))
)
set.seed(1)
out <- ssa(
  initial_state = initial_state,
  reactions = reactions,
  params = params,
  final_time = final_time,
  method = ssa_exact(),
  sim_name = sim_name
)
plot_ssa(out)
