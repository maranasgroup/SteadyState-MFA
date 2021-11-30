function [ mfamodel ] = loadexptdata( mfamodel,filename,expname,append_flag )
%[MFAMODEL] = LOADEXPTDATA(MFAMODEL,FILENAME,EXPNAME,APPEND) reads the
%             experimental data contained within FILENAME and creates a
%             DATA structure to be included in MFAMODEL. the DATA structure
%             contains a list of measured FLUXES, MSDATA, POOL SIZES, and
%             the TRACER scheme for the labeling experiment EXPTNAME. The
%             MSDATA may be corrected for natural abundance of unlabeled
%             atoms based on options specified in MFAMODEL. APPEND allows
%             integration of multiple data sets for improved flux
%             estimation accuracy based on COMPLETE-MFA.
%   Written by Saratram Gopalakrishnan on 09/09/2016 at 102 MRL, Penn State

flux = struct('name',[],'val',[],'std',[],'on',[]);
pool = struct('name',[],'val',[],'std',[],'on',false);
msdata = struct('name',[],'emu',[],'met',[],'mdv',[],'std',[],'mdvreg',[],'on',true);
tracer = struct('met',[],'trc',[]);
trc = struct('name',[],'frac',[],'atmdv',[]);
data = struct('name',expname,'flux',[],'pool',[],'msdata',[],'tracer',[],'on',true);

% Reading expts file

[~,~,raw1] = xlsread(filename,'MSData');
[~,~,raw2] = xlsread(filename,'FluxData');
[~,~,raw3] = xlsread(filename,'TracerData');

% fluxes
flxs = raw2(2:end,:);
for i = 1:length(flxs(:,1))
    f1 = flux;
    f1.name = flxs(i,1);
    f1.val = flxs{i,2};
    f1.std = flxs{i,3};
    f1.on = true;
    data.flux = [data.flux,f1];
end

% pools
data.pool = pool;
if ~mfamodel.options.ss
    [~,~,raw4] = xlsread(filename,'PoolData');
    raw4 = raw4(2:end,:);
    for i = 1:length(raw4(:,1))
        p1 = poo1;
        p1.name = raw4(i,1);
        p1.val = raw4{i,2};
        p1.std = raw4{i,3};
        p1.on = true;
        data.pool = [data.pool,p1];
    end
end

% tracers
if mfamodel.options.sim_na
    defc = [0.9892,0.0108];
else
    defc = [1,0];
end

u = unique(raw3(2:end,1));
for i = 1:length(u)
    tr1 = tracer;
    tr1.met = u(i);
    tx = raw3(ismember(raw3(:,1),u(i)),:);
    f = 0;
    for j = 1:length(tx(:,1))
        tc1 = trc;
        tc1.name = tx(j,2);
        tc1.frac = tx{j,3};
        f = f + tc1.frac;
        nc = tx{j,4};
        tc1.atmdv = zeros(nc,2);
        lpos = tx{j,5};
        if ischar(lpos)
            lpos = str2num(lpos);
        end
        pur = tx{j,6}/100;
        for k = 1:nc
            if ismember(k,lpos)
                tc1.atmdv(k,:) = [(1-pur),pur];
            else
                tc1.atmdv(k,:) = defc;
            end
        end
    end
    tr1.trc = [tr1.trc,tc1];
    if f < 1
        tc1 = trc;
        tc1.name = {['U-13C ',u{i}]};
        tc1.frac = 1-f;
        tc1.atmdv = zeros(nc,2);
        for k = 1:nc
            tc1.atmdv(k,:) = defc;
        end
        tr1.trc = [tr1.trc,tc1];
    end
    data.tracer = [data.tracer,tr1];
end

% MS Data
MDV = raw1(2:end,:);
for i = 1:length(MDV(:,1))
    % Extract data (name, emu, MDV)
    m1 = msdata;
    m1.name = MDV(i,1);
    m1.emu = MDV(i,2);
    [m1.met{1},at] = emu2info(m1.emu{1}); % Get name and MDV from imported data
    mz = str2num(MDV{i,4});
    sz = str2num(MDV{i,5});
    
    % Count atoms from formula
    % For GC-MS fragments, make sure to get the fragment formula
    % For LC-MS, ???
    atms = decode_fragformula(char(MDV{i,3}));
    l = length(at);
    
    % Get atoms that are outside of 
    atms(1) = atms(1) - l;
    
    m1.mdvreg = gen_nat_dil(atms);
    
    % fcor correction
    if mfamodel.options.fcor
        m1.std = sz(1:l+1);
        uc = m1.mdvreg(1:length(mz));
        M = zeros(length(uc));
        M(1,:) = uc;
        M = M';
        for k = 2:length(uc)
            M(k:end,k) = M(k-1:end-1,k-1);
        end
        cr = M\(mz');
        %if length(cr)<(l+1)
        %    cr((length(cr)+1):(l+1)) = 0;
        %end
        cr = cr/sum(cr(1:l+1));
        m1.mdv = [cr(1:l+1)]';
        m1.mdvreg = 1;
    else
        m1.std = sz;
        m1.mdv = mz;
    end
    if mfamodel.options.default_sd
        mx = m1.mdv;
        ms = min(0.5,max(0.005,mx));
        m1.std = 0.003 + ((0.005-0.003)/(0.5-0.005))*(ms-0.005);
    end
    data.msdata = [data.msdata,m1];
end
if append_flag
    mfamodel.data = [mfamodel.data,data];
else
    mfamodel.data = data;
end
mfamodel.level = 3;
end

function [name,atom_ids] = emu2info(emu_id)
% Helper Function to split EMU into metabolites and atom positions

emu = char(emu_id);
pos = strfind(emu,'-');
name = emu(1:(pos(end)-1));
atoms = emu((pos(end)+1):end);
atoms = regexp(atoms,',','split');
atoms = strtrim(atoms);
atom_i = char(atoms);
atom_ids = str2num(atom_i);
if size(atom_ids,1) > 1
    atom_ids = atom_ids';
end
end

function y = gen_nat_dil(atm_f)
% Helper Function to generate labeling distribution based on natural
% abundance

at{1,1} = [0.989184 0.010816];                  %c
at{2,1} = [0.999844 0.000156];                  %h
at{3,1} = [0.99758 0.00038 0.00204];            %o
at{4,1} = [0.996337 0.003663];                  %n
at{6,1} = [0.9222 0.0469 0.0309];               %si
at{5,1} = [0.95039 0.0075 0.0421 0 0.0002];     %s
y = 1;
for i = 1:1:length(atm_f)
    for j = 1:1:atm_f(i)
        y = conv2(y,at{i,1});
    end
end
end
% Helper function to convert fragment formula into atom types
function y = decode_fragformula(string)
atoms = {'c','h','o','n','s','si'};
y = zeros(1,length(atoms));
a = string;
for i = 1:1:length(y)
a = regexp(a,char(atoms(length(y)+1-i)),'split');
if length(a) > 1
    if isempty(char(a(2)))
        a{2} = '1';
    end
    y(length(y) + 1 - i) = str2num(char(a(2)));
end
a = char(a(1));
end
end
