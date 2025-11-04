clear
clc

%{
% Sample data: Concentrations of A and B and observed reaction rate
A_conc = [0.1, 0.2, 0.1, 0.3, 0.2, 0.3, 0.4, 0.5]; % Concentration of A (mol/L)
B_conc = [0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4]; % Concentration of B (mol/L)
rate = [0.15, 0.25, 0.22, 0.35, 0.40, 0.45, 0.50, 0.55]; % Observed reaction rate (mol/L·s)

% Define the Langmuir-Hinshelwood rate equation with element-wise multiplication
langmuir_hinshelwood_eqn = @(params, xdata) (params(1) * params(2) .* xdata(:,1) .* params(3) .* xdata(:,2)) ./ (1 + params(2) .* xdata(:,1) + params(3) .* xdata(:,2)).^2;

% Initial guess for the parameters: [k, K_A, K_B]
initial_guess = [1, 1, 1];

% Combine concentrations A and B into a single matrix for convenience
xdata = [A_conc', B_conc']; % Each row of `xdata` corresponds to [A, B]

% Perform nonlinear regression using 'lsqcurvefit' to minimize the sum of squared errors
options = optimset('Display', 'off'); % Turn off display during optimization
[params_est, resnorm, residual, exitflag, output] = lsqcurvefit(langmuir_hinshelwood_eqn, initial_guess, xdata, rate', [], [], options);

% Extract estimated parameters
k_est = params_est(1);
K_A_est = params_est(2);
K_B_est = params_est(3);

% Display the estimated parameters
fprintf('Estimated Parameters:\n');
fprintf('Rate constant (k) = %.4f\n', k_est);
fprintf('Adsorption constant for A (K_A) = %.4f\n', K_A_est);
fprintf('Adsorption constant for B (K_B) = %.4f\n', K_B_est);

% Plot the observed vs. predicted rates
predicted_rate = langmuir_hinshelwood_eqn(params_est, xdata);

% Check the dimensions of the observed and predicted rates
disp(size(rate));  % Should be [8, 1] or [1, 8]
disp(size(predicted_rate));  % Should be [8, 1] or [1, 8]

figure;
plot3(A_conc, B_conc, rate, 'o', 'MarkerFaceColor', 'b'); % Observed data
hold on;
plot3(A_conc, B_conc, predicted_rate, 'r-', 'LineWidth', 2); % Fitted data
xlabel('[A] (mol/L)');
ylabel('[B] (mol/L)');
zlabel('Rate (mol/L·s)');
legend('Observed Data', 'Fitted Model');
title('Langmuir-Hinshelwood Kinetics Fit');
grid on;
%}

%{
% Sample data: Concentrations of A and B and observed reaction rate
A_conc = [0.1, 0.2, 0.1, 0.3, 0.2, 0.3, 0.4, 0.5]; % Concentration of A (mol/L)
B_conc = [0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4]; % Concentration of B (mol/L)
rate = [0.15, 0.25, 0.22, 0.35, 0.40, 0.45, 0.50, 0.55]; % Observed reaction rate (mol/L·s)

% Define the Langmuir-Hinshelwood-Hougen-Watson (LHHW) rate equation
LHHW_eqn = @(params, xdata) (params(1) * (params(2) * xdata(:,1)) .* (params(3) * xdata(:,2))) ./ ...
                            (1 + params(2) * xdata(:,1) + params(3) * xdata(:,2) + ...
                             params(2) * params(3) * xdata(:,1) .* xdata(:,2)).^2;

% Initial guess for the parameters: [k, K_A, K_B]
initial_guess = [1, 1, 1];

% Combine concentrations A and B into a single matrix for convenience
xdata = [A_conc', B_conc']; % Each row of `xdata` corresponds to [A, B]

% Perform nonlinear regression using 'lsqcurvefit' to minimize the sum of squared errors
options = optimset('Display', 'off'); % Turn off display during optimization
[params_est, resnorm, residual, exitflag, output] = lsqcurvefit(LHHW_eqn, initial_guess, xdata, rate', [], [], options);

% Extract estimated parameters
k_est = params_est(1);
K_A_est = params_est(2);
K_B_est = params_est(3);

% Display the estimated parameters
fprintf('Estimated Parameters:\n');
fprintf('Rate constant (k) = %.4f\n', k_est);
fprintf('Adsorption constant for A (K_A) = %.4f\n', K_A_est);
fprintf('Adsorption constant for B (K_B) = %.4f\n', K_B_est);

% Plot the observed vs. predicted rates
predicted_rate = LHHW_eqn(params_est, xdata);

% Check the dimensions of the observed and predicted rates
disp(size(rate));  % Should be [8, 1] or [1, 8]
disp(size(predicted_rate));  % Should be [8, 1] or [1, 8]

figure;
plot3(A_conc, B_conc, rate, 'o', 'MarkerFaceColor', 'b'); % Observed data
hold on;
plot3(A_conc, B_conc, predicted_rate, 'r-', 'LineWidth', 2); % Fitted data
xlabel('[A] (mol/L)');
ylabel('[B] (mol/L)');
zlabel('Rate (mol/L·s)');
legend('Observed Data', 'Fitted Model');
title('Langmuir-Hinshelwood-Hougen-Watson Kinetics Fit');
grid on;
%}

%{
% Sample data: Concentrations of A and B and observed reaction rate
A_conc = [0.1, 0.2, 0.1, 0.3, 0.2, 0.3, 0.4, 0.5]; % Concentration of A (mol/L)
B_conc = [0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4]; % Concentration of B (mol/L)
rate = [0.15, 0.25, 0.22, 0.35, 0.40, 0.45, 0.50, 0.55]; % Observed reaction rate (mol/L·s)

% Define the Generalized Langmuir-Hinshelwood-Hougen-Watson (GLHHW) rate equation
GLHHW_eqn = @(params, xdata) ...
    (params(1) * prod(params(2:end) .* xdata, 2)) ./ ...
    (1 + sum(params(2:end) .* xdata, 2) + interaction_terms(params(2:end), xdata));

% Interaction term function to compute the sum of pairwise interactions
function int_terms = interaction_terms(K, xdata)
    n = size(xdata, 1); % Number of data points
    m = length(K); % Number of reactants
    
    % Initialize the interaction term
    int_terms = 0;
    
    % Calculate pairwise interaction terms (assuming 2 reactants, extendable to more)
    for i = 1:m
        for j = i+1:m
            % Product of concentrations and adsorption constants
            int_terms = int_terms + K(i) * K(j) * sum(xdata(:,i) .* xdata(:,j)); 
        end
    end
end

% Initial guess for the parameters: [k, K_A, K_B, K_C, ...] (for 2 reactants: [k, K_A, K_B])
initial_guess = [1, 1, 1];  % For two reactants A and B

% Combine concentrations A and B into a single matrix for convenience
xdata = [A_conc', B_conc']; % Each row of `xdata` corresponds to [A, B]

% Perform nonlinear regression using 'lsqcurvefit' to minimize the sum of squared errors
options = optimset('Display', 'off'); % Turn off display during optimization
[params_est, resnorm, residual, exitflag, output] = lsqcurvefit(GLHHW_eqn, initial_guess, xdata, rate', [], [], options);

% Extract estimated parameters
k_est = params_est(1);
K_A_est = params_est(2);
K_B_est = params_est(3);

% Display the estimated parameters
fprintf('Estimated Parameters:\n');
fprintf('Rate constant (k) = %.4f\n', k_est);
fprintf('Adsorption constant for A (K_A) = %.4f\n', K_A_est);
fprintf('Adsorption constant for B (K_B) = %.4f\n', K_B_est);

% Plot the observed vs. predicted rates
predicted_rate = GLHHW_eqn(params_est, xdata);

% Check the dimensions of the observed and predicted rates
disp(size(rate));  % Should be [8, 1] or [1, 8]
disp(size(predicted_rate));  % Should be [8, 1] or [1, 8]

figure;
plot3(A_conc, B_conc, rate, 'o', 'MarkerFaceColor', 'b'); % Observed data
hold on;
plot3(A_conc, B_conc, predicted_rate, 'r-', 'LineWidth', 2); % Fitted data
xlabel('[A] (mol/L)');
ylabel('[B] (mol/L)');
zlabel('Rate (mol/L·s)');
legend('Observed Data', 'Fitted Model');
title('Generalized Langmuir-Hinshelwood-Hougen-Watson Kinetics Fit');
grid on;


% Initial guess for the parameters: [k, K_A] (for 1 reactant: [k, K_A])
initial_guess = [1, 1];  % For one reactant A

% Combine concentrations A into a single matrix for convenience
xdata = A_conc'; % Each row of `xdata` corresponds to [A]

% Perform nonlinear regression using 'lsqcurvefit' to minimize the sum of squared errors
options = optimset('Display', 'off'); % Turn off display during optimization
[params_est, resnorm, residual, exitflag, output] = lsqcurvefit(GLHHW_eqn, initial_guess, xdata, rate', [], [], options);

% Extract estimated parameters
k_est = params_est(1);
K_A_est = params_est(2);

fprintf('Estimated Parameters:\n');
fprintf('Rate constant (k) = %.4f\n', k_est);
fprintf('Adsorption constant for A (K_A) = %.4f\n', K_A_est);

% Predicted rate
predicted_rate = GLHHW_eqn(params_est, xdata);

% Plot
figure;
plot(A_conc, rate, 'bo', 'MarkerFaceColor', 'b');
hold on;
plot(A_conc, predicted_rate, 'r-', 'LineWidth', 2);
xlabel('[A] (mol/L)');
ylabel('Rate (mol/L·s)');
legend('Observed', 'Fitted');
title('GLHHW Fit — Single Reactant');
grid on;
%}

% Sample data: Concentration of A and observed reaction rate
A_conc = [0.1, 0.2, 0.3, 0.4, 0.5]; % Concentration of A (mol/L)
rate = [0.15, 0.25, 0.35, 0.45, 0.50]; % Observed reaction rate (mol/L·s)

% Define simplified GLHHW rate equation for ONE reactant
GLHHW_eqn = @(params, A) (params(1) * params(2) .* A) ./ (1 + params(2) .* A);
% params(1) = k, params(2) = K_A

% Initial guess for parameters [k, K_A]
initial_guess = [1, 1];

% Perform nonlinear regression using lsqcurvefit
options = optimset('Display', 'off');
[params_est, resnorm, residual, exitflag, output] = lsqcurvefit(GLHHW_eqn, initial_guess, A_conc', rate', [], [], options);

% Extract estimated parameters
k_est = params_est(1);
K_A_est = params_est(2);

% Display estimated parameters
fprintf('Estimated Parameters:\n');
fprintf('Rate constant (k) = %.4f\n', k_est);
fprintf('Adsorption constant for A (K_A) = %.4f\n', K_A_est);

% Predicted rates
predicted_rate = GLHHW_eqn(params_est, A_conc');

% Plot observed vs. predicted rates
figure;
plot(A_conc, rate, 'bo', 'MarkerFaceColor', 'b');
hold on;
plot(A_conc, predicted_rate, 'r-', 'LineWidth', 2);
xlabel('[A] (mol/L)');
ylabel('Rate (mol/L·s)');
legend('Observed Data', 'Fitted Model');
title('GLHHW Model Fit (Single Reactant)');
grid on;
