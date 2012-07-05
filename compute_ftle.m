function ftle_ = compute_ftle(max_eigenvalue,timespan)
%FTLE Finite-time Lyapunov exponent

ftle_ = .5*log(max_eigenvalue)/timespan;

end

