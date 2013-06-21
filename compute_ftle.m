function ftle = compute_ftle(max_eigenvalue,timespan)
%FTLE Finite-time Lyapunov exponent

ftle = .5*log(max_eigenvalue)/timespan;


