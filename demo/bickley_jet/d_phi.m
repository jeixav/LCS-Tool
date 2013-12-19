% Forced-damped Duffing oscillator used with aperiodic forcing
function dPhi = d_phi(tau,phi)
dPhi(2,1) = nan;

dPhi(1) = phi(2);
dPhi(2) = -.1*phi(2) - phi(1)^3 + 11*cos(tau);
