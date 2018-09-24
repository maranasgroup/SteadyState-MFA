function [ mfamodel ] = defopt( mfamodel )
% DEFOPT creates a default options set for MFA simulations
opt = struct('ss',true,'sim_na',true,'fcor',true,'default_sd',false,'output_display',true,'dfbase',1e-6,'conf_lvl',0.95,'multistart',1,'reinit',true);
opt.simsens = true; %simulate sensitivities
opt.conf_set = 'minset_main';   % options are: 'all', 'main', 'dilution', 'all_net', 'all_exch', 'minset_main', 'minset_all', 'custom'
opt.conf_custom = [];
opt.conf_step = 5; %expected number of steps required to reach threshold chi2 value
mfamodel.options = opt;

end

