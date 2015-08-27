function fig7_fisherinfo
% Fig. 7 - Fisher information under fluctuations of attentional state


% set input noise level for given threshold
criterion = 0.75;               % threshold criterion (% correct)
dp = 2 * norminv(criterion);    % equivalent d'
thresh = [0.5 3] / 180 * pi;    % threshold based on input noise
sdtheta = thresh / dp;          % SD of input noise

n = 10 : 10 : 500;              % number of neurons
sdalpha = 0.1;                  % variance of spatial gain
mbeta = 0.2;                    % average feature gain
sdbeta = 0.1;                   % variance of feature gain
sdpsi = 10 / 180 * pi;          % variance of attended direction

% tuning parameters
kappa = 2;                      % tuning width
gamma = log(30) - kappa;        % average peak firing rate
amplvar = 1.5;                  % variance of amplitudes across the population

nk = numel(n);
nt = numel(sdtheta);
ns = numel(sdalpha);
J = zeros(nk, ns, nt);

% generate population
for k = 1 : nk
    neurons = n(k);
    phi = (0.5 : neurons) / neurons * 2 * pi - pi;
    f = exp(kappa .* cos(phi) + gamma);
    df = -kappa .* sin(-phi) .* f;
    h = cos(phi);
    lambda = (1 + mbeta * h) .* f;
    dlambda = (1 + mbeta * h) .* df;

    % (1) no attentional fluctuations
    J0 = dlambda * diag(1 ./ lambda) * dlambda';
    J(k, 1, 1) = J0;
    
    % (2) spatial gain fluctuations
    J(k, 2, 1) = J0 - amplvar * sum(dlambda .^ 2) / (sdalpha ^ -2 + sum(f));
    
    % (3) feature gain fluctuations
    J(k, 3, 1) = J0 - amplvar * sum((h .* dlambda) .^ 2) / (sdbeta ^ -2 + sum(h .^ 2 .* f));
    
    % (4) attended feature fluctuations
    J(k, 4, 1) = J0 / (1 + (sdpsi ^ 2 * mbeta ^ 2 / kappa ^ 2) * J0);
    
    % add effect pf input noise
    for j = 1 : nt
        for i = 1 : 4
            J(k, i, j + 1) = J(k, i, 1) / (1 + sdtheta(j) ^ 2 * J(k, i, 1));
        end
    end
end

fig = Figure(7, 'size', [110 50]);
colors = [0 0 0; 0 0.7 0; 0.9 0.7 0.1; 1 0.2 0.2];
set(fig, 'DefaultAxesColorOrder', repmat(colors(1 : 4, :), 2, 1))
Jmax = 5000;
nmax = 500;


% [A] Fisher information in the absence of input noise
subplot(1, 2, 1)
hold on
plot(n, J(:, :, 1))
Jthresh = kappa ^ 2 / (sdpsi ^ 2 * mbeta ^ 2);
plot(xlim, [1 1] * Jthresh, '--k')
text(nmax, Jthresh, sprintf('%.2g', dp ./ sqrt(Jthresh) / pi * 180), 'FontSize', 14)
set(gca, 'xlim', [0 nmax], 'ylim', [0 Jmax])
xlabel('Number of neurons')
ylabel('Fisher information')
legend({'No fluctuations', 'Spatial gain', 'Feature gain', 'Attended feature'})
axis square


% [B] Fisher information with different levels of input noise
subplot(1, 2, 2)
hold on
Jthresh = zeros(2);
plot(n, J(:, 1, 1), ':k')
for i = 1 : nt
    plot(n, J(:, :, i + 1))
    for j = 1 : 2
        Jthresh(i, j) = kappa ^ 2 / (kappa ^ 2 * sdtheta(i) ^ 2 + (j - 1) * sdpsi ^ 2 * mbeta ^ 2);
        plot(xlim, [1 1] * Jthresh(i, j), '--k')
    end
end
set(gca, 'xlim', [0 nmax], 'ylim', [0 Jmax])
xlabel('Number of neurons')
for i = 1 : numel(Jthresh)
    text(nmax, Jthresh(i), sprintf('%.2g', dp ./ sqrt(Jthresh(i)) / pi * 180), 'FontSize', 14)
end
axis square


fig.cleanup();
fig.save([mfilename('fullpath'), '_.eps'])
