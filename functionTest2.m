function [] = functionTest2()
gammaAlpha = (logspace(0,5) - 1) / 1e4; %logspace(a,b):generate 50 values from 10^a to 10^b
gammaAlpha = [- gammaAlpha(end:-1:1), gammaAlpha];

options = optimoptions('fmincon');
options.Algorithm = 'interior-point';
options.Display = 'none';
options.DiffMinChange = 1e-12;
options.DiffMaxChange = 1e-2;
options.TolFun = 1e-8;
options.TolX = 1e-8;
options.UseParallel = false;

Problem.options = options;
Problem.solver = 'fmincon';
Problem.x0 = [0.005, 0.003];
Problem.lb = [1e-6,   1e-6];%lower bounds
Problem.ub = [1-1e-6, 1e-1];%upper bounds
Problem.objective = @(para) sum((piTanH(gammaAlpha) - piPowerlaw(gammaAlpha, para(1), para(2))).^4);

optPara = fmincon(Problem);
optPara = optPara;
disp(optPara);

close all;
figure;
plot(gammaAlpha, piTanH(gammaAlpha));
hold on; 
plot(gammaAlpha, piPowerlaw(gammaAlpha, optPara(1), optPara(2)), '*');
hold on;
plot(gammaAlpha, piPowerlaw(gammaAlpha,0.0645,0.0049))
legend('tanH', 'Powerlaw(opt)','Powerlaw(Bardella Para)');

end

function pi=piTanH(gammaAlpha)
pi0 = 50;
gamma0Dot = 0.01;

pi = pi0 * tanh(gammaAlpha / gamma0Dot);
end

function pi=piPowerlaw(gammaAlpha,p,gamma0dot)
pi0 = 50;
if nargin < 2
    p = 0.05;

end
if nargin < 3
    gamma0Dot = 0.01 / 3;
end
pi =  sign(gammaAlpha).* pi0.*(abs(gammaAlpha) / gamma0dot).^p ;
end