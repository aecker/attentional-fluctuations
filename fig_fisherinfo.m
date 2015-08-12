function fig_fisherinfo
% Fig. XX - Fisher information under fluctuations of attentional state


% set input noise level for given threshold
criterion = 0.75;                   % threshold criterion (% correct)
Jthresh = [0.3 3] / 180 * pi;       % threshold based on input noise
dpthresh = 2 * norminv(criterion);
stdtheta = Jthresh / dpthresh;      % SD of input noise


n = 10 : 10 : 500;
n = 10 .^ (1 : 0.2 : 3);
nk = numel(n);
sigma = 0.1;                        % variance of spatial gain
nu = 0.2;                           % feature gain
tau = 0.1;                          % variance of feature gain
q = 10 / 180 * pi;                  % variance of attended direction

% tuning width
kappa = 2;

% amplitude variability
gamma = log(30) - kappa;            % average peak firing rate
amplvar = 1.5;                     % variance of amplitudes across the population
% mgamma = -0.5 * log(1 + amplvar);   % E[dgamma] s.t. E[exp(dgamma)] = 1
% stdgamma = sqrt(log(1 + amplvar));  % Var[dgamma] s.t. Var[exp(gamma)] = amplvar

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
    J0(1) = dlambda * diag(1 ./ lambda) * dlambda';
    
    % (2) spatial gain fluctuations
    J0(2) = J0(1) - amplvar * sum(dlambda .^ 2) / (sigma ^ -2 + sum(f));
    
    % (3) feature gain fluctuations
    J0(3) = J0(1) - amplvar * sum((h .* dlambda) .^ 2) / (tau ^ -2 + sum(h .^ 2 .* f));
    
    % (4) attended feature fluctuations
    J0(4) = J0(1) / (1 + (q ^ 2 * nu ^ 2 / kappa ^ 2) * J0(1));
        
    for j = 1 : nt
        for i = 1 : 4
            J(k, i, j) = J0(i) / (1 + stdtheta(j) ^ 2 * J0(i));
        end
    end
end

fig = Figure(11, 'size', [50 50]);
colors = [0 0 0; 0 0.7 0; 0.9 0.7 0.1; 1 0.2 0.2];
set(fig, 'DefaultAxesColorOrder', repmat(colors(1 : 4, :), 2, 1))
Jmax = 5000;
% Jmax = 10000;
nmax = 500;
set(gca, 'xlim', [0 nmax], 'ylim', [0 Jmax])
% set(gca, 'xscale', 'log', 'yscale', 'log', 'xlim', n([1 end]), 'ylim', [0 Jmax])
hold on
Jthresh = zeros(2);
for i = 1 : nt
    plot(n, J(:, :, i))
    for j = 1 : 2
        Jthresh(i, j) = kappa ^ 2 / (kappa ^ 2 * stdtheta(i) ^ 2 + (j - 1) * q ^ 2 * nu ^ 2);
        plot(xlim, [1 1] * Jthresh(i, j), '--k')
    end
end
xlabel('Number of neurons')
ylabel('Fisher information')
legend({'No fluctuations', 'Spatial gain', 'Feature gain', 'Attended feature'})
for i = 1 : numel(Jthresh)
    text(nmax, Jthresh(i), sprintf('%.2g', dpthresh ./ sqrt(Jthresh(i)) / pi * 180), 'FontSize', 14)
end

fig.cleanup();
fig.save([mfilename('fullpath'), '_.eps'])
