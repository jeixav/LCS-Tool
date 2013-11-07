%FTLE Calculate Finite-time Lyapunov exponent
function ftle_ = ftle(max_eigenvalue,timespan)

ftle_ = .5*log(max_eigenvalue)/timespan;


