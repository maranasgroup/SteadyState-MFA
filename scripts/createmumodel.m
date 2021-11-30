function [emumodel,emus,mfamodel] = createmumodel(model_filename,atom_mapping_filename,exptdata_filename)
% EMUMODEL = CREATEMUMODEL(F1,F2,F3) creates the EMUMODEL data object
% containing the data required for flux estimation using steady-state
% 13C-MFA. The inputs F1, F2, F3 correspond to the excel files containing
% the model, atom mapping information, and experimental data, respectively.
% These filenames are supplied as a cell object rather than strings. 
% F3 can be an array of cells that contains multiple datasets for simultaneous
% fitting.
% v1.0: all inputs are strings, no support for integrating multiple
% datasets.
% v2.0: All inputs are to be in cell format to allow integration of
% multiple datasets.

% Written by Saratram Gopalakrishnan on 01/28/2018 in 102 MRL Building


model = xls2MFAmodel(model_filename{1});
[~,mfamodel] = includemapping(model,atom_mapping_filename{1});
mfamodel = defopt(mfamodel);
mfamodel = loadexptdata(mfamodel,exptdata_filename{1},'expt1',0);  %0 if fitting corrected data, 1 otherwise
if length(exptdata_filename)>1
    for i = 2:length(exptdata_filename)
        mfamodel = loadexptdata(mfamodel,exptdata_filename{i},['expt',num2str(i)],1);
    end
end
[emumodel,emus] = emutracer(mfamodel);
emumodel.minset = minconfset(emumodel);
end

