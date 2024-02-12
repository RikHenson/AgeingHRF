function [exp_var1,exp_var2] = kernel_variance(DCM)

if ~isfield(DCM,'plot')
    DCM.plot = 0;
end

M  = DCM.M;  % Model
Ep = DCM.Ep; % Parameters

%% Convert the model to a bilinear form
% This helps because there is a closed-form expression for 
% the Volterra kernels of a bilinear system
[M0,M1,L1,L2] = spm_bireduce(M,Ep);

%% Calculate the Volterra kernels
[H0,H1,H2] = spm_kernels(M0,M1,L1,L2,M.N,M.dt);

%% Reconstruct the signal (at the sampling frequency of the inputs)

% Settings
m = 1; % input condition
r = 1; % output region

% Get the Volterra kernels
V0 = H0(r);
V1 = H1(:,r,m);
V2 = H2(:,:,r,m,m);

% Get the input timeseries
u = full(DCM.U.u(:,m));

% Convolve
[y, G1, G2] = reconstruct_signal(u, V1, V2);

%% Plot
if DCM.plot
    x = (0:(length(y)-1)) .* DCM.U.dt;
    figure;
    plot(x,G1,'b');
    hold on
    plot(x,G1+G2,'r');
    legend({'1st order only','1st+2nd order'});
    xlabel('Time (secs)');
end

%% Report total explained variance
ss_var1 = sum(G1.^2);
ss_var2 = sum(G2.^2);
ss_total = ss_var1+ss_var2;

exp_var1 = ss_var1/ss_total*100;
exp_var2 = ss_var2/ss_total*100;

cond_name = DCM.U.name{m};
fprintf('Of the total variance explained by the Volterra kernels (up to the second order) for the %s condition, %2.2f%% was from the 1st kernel and %2.2f%% from the 2nd kernel\n',cond_name,exp_var1,exp_var2);
