function run_fluxconf_func(res,emod)
    %%% Flux confidence interval estimation %%%
    %%% FOR RUN ON CLUSTER, SETUP GUROBI %%%
    run('load_settings.m');

    addpath('/storage/home/hvd5034/Softwares/gurobi9.0.0_linux64/gurobi900/linux64/matlab')
    gurobi_setup()

    tic
    fprintf('Start flux confidence interval estimation\n')
    emod.minset = minconfset(emod);
    [res, impres] = confintestimate(res, emod);

    % Save
    save(strcat(path_expmt, 'fluxconf.mat'), 'res', 'impres');
    toc
end
