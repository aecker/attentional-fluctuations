function fig6_fisherinfo
% Fig. 6 - Fisher information under fluctuations of attentional state


% set input noise level for given threshold
criterion = 0.75;                   % threshold criterion (% correct)
dpcrit = 2 * norminv(criterion);    % equivalent d'
thresh = [0.5 3] / 180 * pi;        % threshold based on input noise
stdtheta = thresh / dpcrit;         % SD of input noise

n = 10 : 10 : 500;                  % number of neurons
sigma = 0.1;                        % variance of spatial gain
nu = 0.2;                           % feature gain
tau = 0.1;                          % variance of feature gain
q = 10 / 180 * pi;                  % variance of attended direction

% tuning parameters
kappa = 2;                          % tuning width
gamma = log(30) - kappa;            % average peak firing rate
amplvar = 1.5;                      % variance of amplitudes across the population

nk = numel(n);
nt = numel(stdtheta);
ns = numel(sigma);
J = zeros(nk, ns, nt);

% generate population
for k = 1 : nk
    neurons = n(k);
    phi = (0.5 : neurons) / neurons * 2 * pi - pi;
    f = exp(kappa .* cos(phi) + gamma);
    df = -kappa .* sin(-phi) .* f;
    h = cos(phi);
    lambda = (1 + nu * h) .* f;
    dlambda = (1 + nu * h) .* df;

    % (1) no attentional fluctuations
    J0 = dlambda * diag(1 ./ lambda) * dlambda';
    J(k, 1, 1) = J0;
    
    % (2) spatial gain fluctuations
    J(k, 2, 1) = J0 - amplvar * sum(dlambda .^ 2) / (sigma ^ -2 + sum(f));
    
    % (3) feature gain fluctuations
    J(k, 3, 1) = J0 - amplvar * sum((h .* dlambda) .^ 2) / (tau ^ -2 + sum(h .^ 2 .* f));
    
    % (4) attended feature fluctuations
    J(k, 4, 1) = J0 / (1 + (q ^ 2 * nu ^ 2 / kappa ^ 2) * J0);
    
    % add effect pf input noise
    for j = 1 : nt
        for i = 1 : 4
            J(k, i, j + 1) = J(k, i, 1) / (1 + stdtheta(j) ^ 2 * J(k, i, 1));
        end
    end
end

fig = Figure(10, 'size', [110 50]);
colors = [0 0 0; 0 0.7 0; 0.9 0.7 0.1; 1 0.2 0.2];
set(fig, 'DefaultAxesColorOrder', repmat(colors(1 : 4, :), 2, 1))
Jmax = 5000;
nmax = 500;

subplot(1, 2, 1)
hold on
plot(n, J(:, :, 1))
Jthresh = kappa ^ 2 / (q ^ 2 * nu ^ 2);
plot(xlim, [1 1] * Jthresh, '--k')
text(nmax, Jthresh, sprintf('%.2g', dpcrit ./ sqrt(Jthresh) / pi * 180), 'FontSize', 14)
set(gca, 'xlim', [0 nmax], 'ylim', [0 Jmax])
xlabel('Number of neurons')
ylabel('Fisher information')
legend({'No fluctuations', 'Spatial gain', 'Feature gain', 'Attended feature'})
axis square

subplot(1, 2, 2)
hold on
Jthresh = zeros(2);
plot(n, J(:, 1, 1), ':k')
for i = 1 : nt
    plot(n, J(:, :, i + 1))
    for j = 1 : 2
        Jthresh(i, j) = kappa ^ 2 / (kappa ^ 2 * stdtheta(i) ^ 2 + (j - 1) * q ^ 2 * nu ^ 2);
        plot(xlim, [1 1] * Jthresh(i, j), '--k')
    end
end
set(gca, 'xlim', [0 nmax], 'ylim', [0 Jmax])
xlabel('Number of neurons')
for i = 1 : numel(Jthresh)
    text(nmax, Jthresh(i), sprintf('%.2g', dpcrit ./ sqrt(Jthresh(i)) / pi * 180), 'FontSize', 14)
end
axis square

fig.cleanup();
fig.save([mfilename('fullpath'), '_.eps'])
