%{
A = importdata('data_Ct.txt', '\t', 1);
texp = A.data(:, 1);
yexp = A.data(:, 2);

fun = @(r) objectivefunction (r, texp, yexp);
m0 = [1, 1];
m = lsqnonlin(fun, m0);

plot(texp, yexp, 'ko', texp, m0(1) * exp(-m0(2) * texp), 'b-', texp, m(1) * exp(-m(2) * texp), 'b-');
legend('Data','Best fit');
xlabel('t');
ylabel('y');
%}


A = importdata('data_CTt.txt', '\t', 1);
tmpexp = A.data(:, 1);
texp = A.data(:, 2);
yexp = A.data(:, 3);

fun = @(r) objectivefunction (r, tmpexp, texp, yexp);
parameters_0 = [1, 1, 1];
options = optimoptions(@lsqnonlin);
options.Algorithm = 'levenberg-marquardt';
parameters = lsqnonlin(fun, parameters_0, [], [], options);

residuals = yexp' + objectivefunction (parameters, tmpexp, texp, yexp);
ssq = sqrt(sum((yexp-residuals').^2));

plot(texp, yexp, 'ko', texp, yexp' + objectivefunction (parameters, tmpexp, texp, yexp), 'b-');


coeff = 0.5;
[t,y] = ode45(@(t,y) odesrhs(t, y, coeff), [0 20], [0.097]);
plot(t,y(:,1),'-o');

