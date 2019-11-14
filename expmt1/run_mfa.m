%%% BUILD MODEL %%%
tic
fprintf('Start building model\n')
% Load stoichiometric model
model=xls2MFAmodel('./model/example_GEM.xlsx');

% Load atom mapping model
[model,mfamodel]=includemapping(model,'./model/example_AMM.xlsx');

% Create default optimization options
[mfamodel]=defopt_expmt1(mfamodel);

% Load experimental data
[mfamodel]=loadexptdata(mfamodel,'./expmt1/expmt1_data.xlsx','expmt1',false);

% Run EMU
[emod,emus]=emutracer(mfamodel);

% Save
save('./expmt1/models.mat', 'model', 'mfamodel', 'emus', 'emod');
toc
%%% Non-linear minimization of SSR %%%
%[res, foptCell] = flxestimate_fast(emod)

tic
fprintf('Start non-linear optimization\n')
[res, foptCell, residualCell] = flxestimate_proper(emod, 100, 0);

% Save
save('./expmt1/res.mat', 'model', 'mfamodel', 'emus', 'emod',...
    'res', 'foptCell', 'residualCell');
toc

%%% Flux confidence interval estimation %%%
tic
fprintf('Start flux confidence interval estimation\n')
emod.minset = minconfset(emod);
[res, impres] = confintestimate(res, emod);

% Save
save('./res_fluxconf.mat', 'model', 'mfamodel', 'emus', 'emod',...
    'res', 'foptCell', 'residualCell', 'impres');
toc