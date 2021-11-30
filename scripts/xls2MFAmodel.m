function [ model ] = xls2MFAmodel( filename )
%   [model] = xls2MFAmodel(filename) parses the excel file
%   containing details about reactions and metabolites to generate a model
%   suitable for 13C-MFA. unlike a typical metabolic model, an MFA model
%   contains unbalanced metabolites representing sources and sinks
%   MODEL contains the following fields:
%   rid:        reaction id
%   RxnName:    reaction name
%   RxnFormula: the actual reaction
%   metprops:   metabolites structure containing the following fields:
%               metid:      metabolite as in the model
%               metname:    the name of each metabolite
%               metformula: chemical formula of the metabolite
%               ubflag:     unbalanced internal metabolite
%               src:        indicates whether a metabolite is a substrate
%               snk:        indicates whether a metabolite is a sink
%               symmap:     symmetry mapping for metabolite. Default is
%                           blank
%   rxnprops:   reactions structure containing the following fields:
%               rev:        reversibility flag
%               macro:      whether the reaction describes macromolecule
%                           synthesis
%               subsystem:  pathway of which this reaction is a part
%               lb:         lower bound (0 if irreversible, -1E8 if
%                           reversible
%               ub:         upper bound (1E8 by default)
%               baseflx:    properties structure describing reactants,
%                           products and stoichiometric coefficients.
%   S:          the matrix of atoichiometric coefficients
%   lb:         lower bound
%   ub:         upper bound
%   level:      1 = FBA only
%               2 = mapping model incorporated
%               3 = experimental data (flux, pool size and MS
%                   measurements) included
%               4 = parametrized and ready for MFA
%
%   Written By Saratram Gopalakrishnan on 09/03/2016 at 102 MRL, Penn State

[~, ~, raw] = xlsread(filename,'Reactions');
[~,~,raw2] = xlsread(filename,'Metabolites');
model = struct('rid',[],'RxnName',[],'RxnFormula',[],'metprop',[],'rxnprop',[],...
        'S',[]);
    
%basic fills
model.rid = raw(2:end,1);
model.RxnName = raw(2:end,2);
model.RxnFormula = raw(2:end,3);

%metabolites section
lmet = size(raw2,1);
[metprop(1:lmet-1,1)] = deal(struct('metid',[],'metname',[],...
                        'metformula',[],'ubflag',false,...
                        'src',false,'snk',false,'symmap',[],'nC',[]));
for i = 2:lmet
    metprop(i-1,1).metid = raw2(i,1);
    metprop(i-1,1).metname = raw2(i,2);
    metprop(i-1,1).metformula = raw2(i,3);
    metprop(i-1,1).nC = 0;
end
model.metprop = metprop;

%reactions section
lrxn = size(raw,1);
[rxnprop(1:lrxn-1,1)] = deal(struct('rev',false,'lb',0,'ub',10000000,...
    'baseflx',struct('reactant',[],'rstoic',[],'product',[],'pstoic',[],...
    'rmap',[],'pmap',[])));
for i = 2:lrxn

    rxnprop(i-1,1).rev = logical(cell2mat(raw(i,6)));
    rxnprop(i-1,1).macro = logical(cell2mat(raw(i,7)));
    rxnprop(i-1,1).subsystem = (raw(i,8));
    rxnprop(i-1,1).lb = cell2mat(raw(i,4));
    rxnprop(i-1,1).ub = cell2mat(raw(i,5));
    rxnprop(i-1,1).baseflx = parserxn(raw{i,3},rxnprop(i-1,1).baseflx);
end
model.rxnprop = rxnprop;

%S-matrix
model.S = zeros(lmet-1,lrxn-1);
for i = 1:1:size(model.S,2)

    b = rxnprop(i,1).baseflx;
    for j = 1:1:length(b.reactant)
        p = find(ismember([metprop.metid]',b.reactant(j)));
        model.S(p,i) = model.S(p,i) - b.rstoic(j);
    end
    for j = 1:1:length(b.product)
        p = find(ismember([metprop.metid]',b.product(j)));
        model.S(p,i) = model.S(p,i) + b.pstoic(j);
    end
end

% identifying sources sinks and unbalanced metabolites
rev = [rxnprop.rev];
for i = 1:1:size(model.S,1)
    f = find(model.S(i,:));
    if length(f) == 1 && ~rev(f)
        if model.S(i,f) < 0
            model.metprop(i).src = true;
        else 
            model.metprop(i).snk = true;
        end
    elseif length(f) == 1 && rev(f)
        model.metprop(i).ubflag = true;
    end
    k1 = find(model.S(i,:)>0);
    k2 = find(model.S(i,:)<0);
    if length(k1) == length(f) && ~any(rev(f)) && length(f) > 1
        model.metprop(i).ubflag = true;
    elseif length(k2) == length(f) && ~any(rev(f)) && length(f) > 1
        model.metprop(i).ubflag = true;
    end
end

% model capabilities and statistics:
model.level = 1;
disp('FBA model constructed')

end

function flxstruct = parserxn(rxnstring,flxstruct)
% helper function for parsing a reaction string and extracting species and
% stoichiometric coefficients information

% identify common compartment if it exists
if rxnstring(1) == '['
    pos1 = strfind(rxnstring,':');
    pos1 = pos1(1);
    comp = strtrim(rxnstring(1:pos1-1));
    rxnstring(1:pos1) = [];
else
    comp = '';
end

%separate reactants and products side
k = strfind(rxnstring,'>');
pdt = strtrim(rxnstring(k+1:end));
reac = rxnstring(1:k-1);
if isempty(strfind(reac,'<'))
    while(reac(end) == '-')
        reac(end) = [];
    end
    reac = strtrim(reac);
else
    k = strfind(reac,'<');
    reac(k:end) = [];
    reac = strtrim(reac);
end


%reactants side
reac = regexp(reac,'+','split');
reac = strtrim(reac);
for i = 1:length(reac)
    % separate met name and stoich coeff
    k1 = strsplit(reac{i},{' ';'*'});
    if length(k1) == 1
        stx = 1;
        flag = false;
    else
        stx = k1{1};
        k1(1) = [];
        flag = true;
    end
    spectemp = {[k1{1},comp]};
    %handling stoich coeff
    if flag
        if isnan(str2double(stx(1)))
            stx(1) = [];
        end
        if isnan(str2double(stx(end)))
            stx(end) = [];
        end
        stx = str2double(stx);
    end
    if ~any(mod(stx,1))
        tstoic = ones(stx,1);
        [rcts(1:stx,1)] = deal(spectemp);
        rmap = true(stx,1);
        if stx == 0
            tstoic = 0;
            rcts = spectemp;
            rmap = true;
        end
    else
        tstoic = stx;
        rcts = spectemp;
        rmap = false;
    end
    flxstruct.reactant = [flxstruct.reactant;rcts];
    flxstruct.rstoic = [flxstruct.rstoic;tstoic];
    flxstruct.rmap = [flxstruct.rmap;rmap];
    rcts = {};
end

%products side
reac = regexp(pdt,'+','split');
reac = strtrim(reac);
for i = 1:length(reac)
    % separate met name and stoich coeff
    k1 = strsplit(reac{i},{' ';'*'});
    if length(k1) == 1
        stx = 1;
        flag = false;
    else
        stx = k1{1};
        k1(1) = [];
        flag = true;
    end
    spectemp = {[k1{1},comp]};
    %handling stoich coeff
    if flag
        if isnan(str2double(stx(1)))
            stx(1) = [];
        end
        if isnan(str2double(stx(end)))
            stx(end) = [];
        end
        stx = str2double(stx);
    end
    if ~any(mod(stx,1))
        tstoic = ones(stx,1);
        [rcts(1:stx,1)] = deal(spectemp);
        pmap = true(stx,1);
        if stx == 0
            tstoic = 0;
            rcts = spectemp;
            pmap = true;
        end
    else
        tstoic = stx;
        rcts = spectemp;
        pmap = false;
    end
    flxstruct.product = [flxstruct.product;rcts];
    flxstruct.pstoic = [flxstruct.pstoic;tstoic];
    flxstruct.pmap = [flxstruct.pmap;pmap];
    rcts = {};
end    


end
