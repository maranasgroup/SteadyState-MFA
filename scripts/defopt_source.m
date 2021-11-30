function [ mfamodel ] = defopt( mfamodel )
%%% DEFOPT: creates a default options set for MFA simulations %%%
% struct
%   ss: steady state
%   sim_na: simulate natural abundances
%       Accounting for naturally occuring labeled carbon atoms...
%       in the metabolite formula (e.g., pyruvate, accounting for...
%       the fraction of naturally occuring 13C of 3 carbons).
%   fcor: fragment correction
%       Accounting for naturally occuring labeled atoms...
%       that are not parts of traced atoms (i.e., for 13C...
%       MFA, the traced atoms are carbons in metabolites). Thus, atoms...
%       such as carbon and other atoms due to derivization will be ...
%       accounted for (in mass shift).
%   default_sd: flag for using default standard deviation for labeling data
%       of 0.05
%   output_display: output display
%   dfbase: damping factor for step-length for non-linear optimization
%       (need to confirm?)
%   conf_lvl: confidence level of goodness-of-fit (chi-squared test)
%   multistart: number of times to perform non-linear minimization of SSR
%
% opt.conf_set: set of reactions to run flux range estimation subject to
%   confidence interval
%       all? main? dilution? all_net? all_exch? minset_main? minset_all?
%       custom?
% opt.conf_custom: ?
% opt.conf_step: expected number of steps required to reach threshold chi2 value
opt = struct('ss',true,...
            'sim_na',false,... %works
            'fcor',false,... %error: index exceeds the number of array elements (1)
            'default_sd',false,...
            'output_display',true,...
            'dfbase',1e-6,...
            'conf_lvl',0.95,...
            'multistart',10);
opt.conf_set = 'minset_main';   % options are: 'all', 'main', 'dilution', 'all_net', 'all_exch', 'minset_main', 'minset_all', 'custom'
opt.conf_custom = [];
opt.conf_step = 10;
mfamodel.options = opt;

end

