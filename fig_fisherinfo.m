function fig_fisherinfo
% Fig. XX - Fisher information under fluctuations of attentional state


% set input noise level for given threshold
criterion = 0.75;                   % threshold criterion (% correct)
Jthresh = [0.3 3] / 180 * pi;       % threshold based on input noise
dpthresh = 2 * norminv(criterion);
stdtheta = Jthresh / dpthresh;      % SD of input noise

nk = 7;
k0 = 2;
runs = 100;
sigma = [0 0.1 0 0];                % variance of spatial gain
nu = 0.2;
tau = [0 0 0.1 0];                  % variance of feature gain
q = [0 0 0 10] / 180 * pi;          % variance of attended direction
mgamma = 2.3;
stdgamma = 1.35;

nt = numel(stdtheta);
ns = numel(sigma);
J = zeros(nk, ns, nt, runs);

for r = 1 : runs
    % generate population
    neurons = 2 ^ (k0 + nk);
    phi = (0.5 : neurons) / neurons * 2 * pi - pi;
    kappa = 2;
    gamma = zscore(randn(1, neurons)) * stdgamma + mgamma - kappa;
    f = exp(kappa .* cos(phi) + gamma);
    df = -kappa .* sin(-phi) .* f;
    for l = 1 : nt
        for i = 1 : ns
            for k = 1 : nk
                Jk = zeros(1, nk - k + 1);
                for j = 1 : nk - k + 1
                    skip = 2 ^ (nk - k);
                    ndx = j : skip : neurons;
                    h = cos(phi(ndx));
                    lambda = (1 + nu * h) .* f(ndx);
                    dlambda = (1 + nu * h) .* df(ndx);
                    Ci = r1upd(diag(1 ./ lambda), stdtheta(l) * dlambda);
                    if sigma(i) > 0
                        Ci = r1upd(Ci, sigma(i) * lambda);
                    end
                    if tau(i) > 0
                        Ci = r1upd(Ci, tau(i) * h .* lambda);
                    end
                    if q(i) > 0
                        dh = -sin(phi(ndx));
                        Ci = r1upd(Ci, q(i) * nu * dh .* lambda);
                    end
                    Jk(j) = dlambda * Ci * dlambda';
                end
                J(k, i, l, r) = mean(Jk);
            end
        end
    end
end


fig = Figure(10, 'size', [50 50]);
colors = [0 0 0; 0 0.7 0; 0.9 0.7 0.1; 1 0.2 0.2];
set(fig, 'DefaultAxesColorOrder', repmat(colors(1 : 4, :), 2, 1))
Jmax = 5000;
n = 500;
set(gca, 'xlim', [0 n], 'ylim', [0 Jmax])
hold on
Jthresh = zeros(2);
for i = 1 : nt
    plot(2 .^ (k0 + (1 : nk)), mean(J(:, :, i, :), 4))
    for j = 1 : 2
        Jthresh(i, j) = kappa ^ 2 / (kappa ^ 2 * stdtheta(i) ^ 2 + q(end + 1 - j) ^ 2 * nu ^ 2);
        plot(xlim, [1 1] * Jthresh(i, j), '--k')
    end
end
xlabel('Number of neurons')
ylabel('Fisher information')
legend({'No fluctuations', 'Spatial gain', 'Feature gain', 'Attended feature'})
for i = 1 : numel(Jthresh)
    text(n, Jthresh(i), sprintf('%.2g', dpthresh ./ sqrt(Jthresh(i)) / pi * 180), 'FontSize', 14)
end

fig.cleanup();
fig.save([mfilename('fullpath'), '_.eps'])



function C = r1upd(Ai, b)
% Matrix inversion by rank-one update.
%   C = inv(A + b*b'), where Ai = inv(A) is known and b is a vector.

Aib = Ai * b(:);
C = Ai - (Aib * Aib') / (1 + b(:)' * Aib);
