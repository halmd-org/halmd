% Reference implementation to compute force, energy, ... for power_law_with_core module


this_is_a_script = 1;

function [r_, fval, en_pot, hvir] = plwc(cutoff, core, epsilon, sigma, n, r)
% -- function [r_, fval, en_pot, hvir] = plwc(cutoff, core, epsilon, sigma, n, r)
% @param cutoff cutoff radius [sigma]
% @param core core radius [sigma]
% @param epsilon interaction strength [MD units]
% @param sigma length scale [MD units]
% @param n power index of interaction [pos. integer]
% @param r distance between particles [MD units]

    % Adapt precision
    [cutoff, core, epsilon, sigma, n, r] = single(cutoff, core, epsilon, sigma, n, r);
    % Calculate the output
    r_ = r;
    pot_cutoff = epsilon .* (1.0 ./ (cutoff - core)).^n;
    en_pot = epsilon .* (sigma ./ (r - sigma.*core)).^n - pot_cutoff;
    fval = n .* en_pot ./ (r .* (r - sigma.*core));
    hvir = n .* en_pot .* r .* (n .* r + sigma.*core) ./ (r - sigma.*core).^2;
end


% output 13 digits in scientific notation
format("long");

% AA interaction
r = [0.5; 1; 3; 5; 10];
[R_, fval, en_pot, hvir] = plwc(5, 0.4, 1, 1, 6, r);
out = [R_, fval, en_pot, hvir]

% AB interaction
r = [1.5; 3; 5; 9.3; 15];
[R_, fval, en_pot, hvir] = plwc(5, 0.5, 0.5, 2, 12, r);
out = [R_, fval, en_pot, hvir]

% BB interaction
r = [1.75; 4; 7; 12.2; 19];
[R_, fval, en_pot, hvir] = plwc(5, 0.6, 0.25, 2.2, 7, r);
out = [R_, fval, en_pot, hvir]
