% Seismic design codes in Chile allow structural desginers to make only 10 Non Linear Time History analyses and from those analyses one is the maximum number of records. 
% Probabilidad de que de 10 "trials" se obtenga 1 éxito -> es como binomial pero con "p" variable para cada registro sísmico (se multiplican en vez de usar exponente)

clear variables
close all
clc

n_trials = 10; % number of trials
p_range_Values = [0.1, 0.2, 0.3, 0.4, 0.35, 0.23]; % this must be the hard thing to predict

% Plot (revisar esta parte )
Trial_legend = cell(n_trials, 1); % Initialize legend cell array
figure
hold on
for n = 1:n_trials
    k_values = (0:n).'; % Initialize number of successes
    y = zeros(n,1); % Initialize probability of success
    for k = 1:length(k_values) % Use n_trials as upper limit
        mult = 1; % Initialize probability of success for each trial
        for i = 1:k_values(k)
            random_p = p_range_Values(randi(length(p_range_Values)));
            random_q = 1 - random_p;
            mult = mult*random_p*random_q; % Update probability of success
        end
        y(k) = nchoosek(n, k)*mult; % Calculate binomial coefficient
    end
    plot(k_values, y)
    Trial_legend{n} = ['Trial ' num2str(n)]; % Create legend for each trial
end
hold off
xlabel('Number of successes (k)')
ylabel('Probability')
title('Probability of different number of collapses in n trials')
legend(Trial_legend)