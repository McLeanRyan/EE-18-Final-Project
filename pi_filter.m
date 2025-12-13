%%Define constants to use
%Frequency either 13.6MHz or 25.45GHz
f = 13.6*10^6;
w = 2 * pi * f;
%Input Impedance
Z_0 = 50;
%Load Impedance, real and imaginary
Z_load =7 - 22.5j;

%impedance calculations:
%Ar @3W = 7 - j22.5
%02 @ 18W = 6.59 - j22.62

%Set fixed inductance value  
L_fixed = 0.4 * 1e-6;
L_fixed_uH = L_fixed * 1e6; % Defined this so the final print statement works
X_L_fixed = w * L_fixed;
fprintf('Fixed Inductor: %.8f H\n', L_fixed);
fprintf('Target Load: %s Ohms\n', num2str(Z_load));
%Convert load to parallel equiv.
Y_load = 1 / Z_load;
G_load = real(Y_load);
B_load = imag(Y_load);
R_p = 1 / G_load;      % Parallel Resistance of Load
X_p_intrinsic = -1/B_load; % Parallel Reactance of Load
%Virtual resistor for pi-match
% Equation: X_L = Rv * ( sqrt(Z0/Rv - 1) + sqrt(Rp/Rv - 1) )
%Solved w fzero
% FIXED: Changed Z0 to Z_0 in the equation handle below
match_eqn = @(Rv) (Rv * (sqrt(Z_0/Rv - 1) + sqrt(R_p/Rv - 1))) - X_L_fixed;
% The Virtual Resistance must be lower than both Z0 and Rp
%0.01 is subtracted to make the program stable / not accidentally root a
%negative
max_Rv = min(Z_0, R_p) - 0.01;
try
    Rv = fzero(match_eqn, [0.1, max_Rv]);
catch
    error('Inductor is missized.');
end
Q1 = sqrt(Z_0/Rv - 1);  % Source Side Q
Q2 = sqrt(R_p/Rv - 1); % Load Side Q
fprintf('Virtual Resistance Rv = %.2f Ohms\n', Rv);
%capacitor calculations
%C1 cancels Q1 reactance
% FIXED: Changed Z0 to Z_0 below
X_C1 = -Z_0 / Q1;
C_tune = -1 / (w * X_C1);
%C2 
%C2 should provide reactance for Q2 and cancel load / plasma reactance
B_target_total = Q2 / R_p;         % Susceptance needed for the match
B_intrinsic = B_load;              % Susceptance provided by plasma
B_external_needed = B_target_total - B_intrinsic;
% Convert B to C
% If B is positive, we need a Capacitor. (B = wC)
C_load = B_external_needed / w;
%% 6. OUTPUT RESULTS
fprintf('\nREQUIRED COMPONENTS:\n');
fprintf('--------------------\n');
fprintf('1. Fixed Inductor (Series):  %.2f uH\n', L_fixed_uH);
fprintf('2. Tune Capacitor (Source):  %.2f pF\n', C_tune * 1e12);
fprintf('3. Load Capacitor (Plasma):  %.2f pF\n', C_load * 1e12);