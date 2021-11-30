function run_minSSR_func()
    %%% BUILD MODEL %%%
    run('load_settings.m');

    tic
    fprintf('Start building model\n')
    % Load stoichiometric model
    model=xls2MFAmodel(gem_file);

    % Load atom mapping model
    [model,mfamodel]=includemapping(model, amm_file);

    % Create default optimization options
    [mfamodel]=defopt(mfamodel);

    % Load experimental data
    if size(expmt_files, 2) == 1
        fpath = strcat(path_expmt, expmt_files{1});
        fid = expmt_ids{1};
        [mfamodel]=loadexptdata(mfamodel, fpath, fid, false);
    else
        for i = 1:size(expmt_files, 2)
            fpath = strcat(path_expmt, expmt_files{i});
            fid = expmt_ids{i};
            [mfamodel]=loadexptdata(mfamodel, fpath, fid, true);
        end
    end

    % Run EMU
    [emod,emus]=emutracer(mfamodel);

    % Save
    save(strcat(path_expmt, 'models.mat'), 'model', 'mfamodel', 'emus', 'emod');
    toc
    %%% Non-linear minimization of SSR %%%
    %[res, foptCell] = flxestimate_fast(emod, randscale)

    tic
    fprintf('Start non-linear optimization\n')
    [res, foptCell, residualCell] = flxestimate_proper(emod, repeat, randscale);

    % Save
    save(strcat(path_expmt, 'res.mat'), 'model', 'mfamodel', 'emus', 'emod',...
        'res', 'foptCell', 'residualCell');
    toc
    
end
