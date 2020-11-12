%% Set input parameters
m = 3;
delta = 0.15;
f_s = 260;

%% Get band edges
q = floor(0.1*m);
r = m - 10*q;
f_sl = 25 + 1.9*q + 4.1*r;
f_sh = f_sl + 20;
f_pl = f_sl - 4;
f_ph = f_sh + 4;

%% Normalize and convert to low-pass analog
f = [f_pl f_sl f_sh f_ph];
w = 2*pi*f/f_s;
sigma = tan(w/2);
sigma_0 = sqrt(sigma(1)*sigma(4));
B = sigma(4) - sigma(1);
sigma_l = (B*sigma)./(sigma_0^2 - sigma.^2);
sigma_lp = 1;
sigma_ls = min(abs([sigma_l(2), sigma_l(3)]));

%% Analog Lowpass Transfer function
D1 = 1/(1-delta)^2 - 1;
epsilon = sqrt(D1);
D2 = 1/delta^2 - 1;
N = ceil(log(sqrt(D2/D1))/log(sigma_ls/sigma_lp));

% since cutoff frequency is between the two expressions, 
% and average is a good number in between to set it to
sigma_c = 0.5*((sigma_lp/(D1^(1/(2*N)))) + (sigma_ls/(D2^(1/(2*N)))));

syms s;
poles = double(solve(1+(s/(1i*sigma_c))^(2*N)==0, s));
% get the poles that lie in the left-half of the plane
admissible_poles = poles(real(poles) < 0);
% construct a polynomial using these poles
denom = poly(admissible_poles);
num = denom(N+1);

syms s z;
H_analog_lpf(s) = poly2sym(num, s)/poly2sym(denom, s);

%% Bandpass and Discrete Transfer Function
H_analog_bpf(s) = H_analog_lpf((s*B)/(s^2 + sigma_0^2));
H_discrete_bpf = H_analog_bpf((z-1)/(z+1));

% get coefficients of analog bpf
[num_analog, denom_analog] = numden(H_analog_bpf);
num_analog = sym2poly(expand(num_analog));
denom_analog = sym2poly(expand(denom_analog));
k = denom_analog(1);
num_analog = num_analog/k;
denom_analog = denom_analog/k;

% get coefficients of discrete bpf
[num_discrete, denom_discrete] = numden(H_discrete_bpf);
num_discrete = sym2poly(expand(num_discrete));
denom_discrete = sym2poly(expand(denom_discrete));
k = denom_discrete(1);
num_discrete = num_discrete/k;
denom_discrete = denom_discrete/k;
[H,f] = freqz(num_discrete, denom_discrete, 1000, 260e3);

%% Required Plots
fvtool(num_discrete, denom_discrete, 'Analysis', 'freq')
figure
hold on
plot(f,abs(H), "DisplayName", "Magnitude Response")
title("Magnitude Response of Bandstop IIR Filter")
xlabel("Frequency (Hz)")
ylabel("Magnitude")
plot(f, ones(size(f))*0.15, '--', "DisplayName", "0.15");
plot(f, ones(size(f))*0.85, '--', "DisplayName", "0.85");
plot(f, ones(size(f))*1.15, '--', "DisplayName", "1.15");
plot([1 1]*f_sl*1e3, [0 1.15], '--', "DisplayName", "f_{s1}");
plot([1 1]*f_pl*1e3, [0 1.15], '--', "DisplayName", "f_{p1}");
plot([1 1]*f_ph*1e3, [0 1.15], '--', "DisplayName", "f_{p2}");
plot([1 1]*f_sh*1e3, [0 1.15], '--', "DisplayName", "f_{s2}");
% plot(f,angle(H), "DisplayName", "Phase Response")
grid
legend
hold off