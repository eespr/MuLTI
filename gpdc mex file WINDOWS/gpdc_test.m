% gpdc -R 5 -s frequency -n 150 -min 1.0 -max 150.0 test.model | figue -c

T  = [7.5, 25, 0];
Vp = [500, 1350, 2000];
Vs = [200, 210, 1000];
d  = [1700, 1900, 2500];
freq = linspace(1,150,150);
out = gpdc(T, Vp, Vs, d, 'fV', freq);

fig00 = figure;
plot(out(:, 1), out(:, 2:end));

fig01 = figure;
plot(out(:, 1), rdivide(1, out(:, 2:end)));

% same thing, using nSamples, minRange and maxRange:
%   out = gpdc(T, Vp, Vs, d, 'minRange', 1.0, 'maxRange', 150.0, 'nSamples', 150);

% using fV (frequency vector):
%   out = gpdc(T, Vp, Vs, d, 'fV', [10, 20, 30]);
% should be equivalent to:
%   gpdc -R 5 -s frequency -n 3 -min 10.0 -max 30.0 test.model

