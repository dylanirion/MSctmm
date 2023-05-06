option for local updates instead of global
Update example
Write better parameter checks (Hmat dims, tau_pos >= tau_vel, etc)
Clean everything up!
sim function doesn't know that OUF states should share mu
covariates for each x,t
 - gamma_j (transitions out of a state) = sum gamma_ij where i != j. (diagonal of rate matrix)
 - gamma_j is globally bounded above
 - kappa = max(gamma_j)
 - potential switch has probability of being an actual switch by gamma_j/kappa
 - simulate state sequence - generate potential switches T1, T2, ... with rate kappa Poisson process
 - for each k, generate location x(T_k) by forward simulation, then decide if it is a true switch with probability gamma_J/kappa
 - if true switch, pick state with probability gamma_J_j/gamma_J
 - retain all x(Tk) as part of the trajectory, but only those true switches will be switches
 gamma_i = alpha/(1+exp(-alpha*(tmod1-t_alpha))) where t_alpha is most likely day of year, alpha is rate parameter, and t is current day of the year
 normal priors of log(alpha) truncated above at log(kappa)
 uniform prior for t_alpha on (0,365.25) expressed in days?
 fix kappa to all for certain number of transitions

Semi Markov? transitions from hazard rate of distribution other than exponential

RJMCMC?
abstract out "model" arg (accept function, same for priors)
