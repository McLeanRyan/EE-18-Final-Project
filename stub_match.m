clear
% Define Set Impedences of the System
ZR = 50;
ZL_o2 = 6.59-22.62j;
ZL_Ar = 7-22.5j;
Z0 = 50;

% Admittances
YR = 1 / ZR;
YL_o2 = 1 / ZL_o2;
YL_Ar = 1 / ZL_Ar;

ZL = ZL_o2;
YL = 1/(ZL/Z0);
Y0 = 1 / Z0;

% Define Frequency Variables
freq1 = 13600000;
freq2 = 2450000000;

lambda1 = (3*10^8)/freq1;
lambda2 = (3*10^8)/freq2;

lambda1 = lambda2;

% Define Power Variables
PO2 = 18;
PAr = 3;

% Define Stub Variables
l1 = 0;
l2 = 0;
d1 = lambda1; %Adjustable
d2 = 0.375*lambda1; %Adjustable

% Calculate Stub Lengths
Y_d1 = (YL + 1j*tan(2*pi*d1)) / (1 + 1j*YL*tan(2*pi*d1));

Y_BB = Y_d1 + Y0;

Y_d2 = (Y_BB + 1j*tan(2*pi*d2)) / (1+1j*Y_BB*tan(2*pi*d2));

Y_AA = Y_d2 + Y0;

Y_AA_re = real(Y_AA);

m = tan(2*pi*d2);

B_s11 = (1-m*imag(Y_d1)-sqrt(real(Y_d1)*((1+m^2)-real(Y_d1)*m^2))) / m;
B_s12 = (1-m*imag(Y_d1)+sqrt(real(Y_d1)*((1+m^2)-real(Y_d1)*m^2))) / m;

B_s21 = (1/m)*(1-sqrt((1+m^2-real(Y_d1)*m^2)/real(Y_d1)));
B_s22 = (1/m)*(1+sqrt((1+m^2-real(Y_d1)*m^2)/real(Y_d1)));

b_s11 = -1 / B_s11;
b_s12 = -1 / B_s12;
n_s21 = -1 / B_s21;
n_s22 = -1 / B_s22;

l_11 = (lambda1 / (2*pi)) * atan(b_s11);
l_12 = (lambda1 / (2*pi)) * atan(b_s12);
l_21 = (lambda1 / (2*pi)) * atan(n_s21);
l_22 = (lambda1 / (2*pi)) * atan(n_s22);

if l_11<0
    l_11 = l_11 + lambda1/2;
end
if l_12<0
    l_12 = l_12 + lambda1/2;
end
if l_21<0
    l_21 = l_21 + lambda1/2;
end
if l_22<0
    l_22 = l_22 + lambda1/2;
end


