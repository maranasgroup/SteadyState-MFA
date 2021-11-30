function [ emod,emus ] = emutracer( mfamodel )
% EMU Decomposition function
%emumodel = mfamodel;

data = mfamodel.data;
nxp = length(data);

%identifying list of metabolites involved in EMU balances
M0 = false(length([mfamodel.mapnode.met]),1);
for i = 1:nxp
    met = [data(i).msdata.met]';  % collecting list of metabolites with fragments to be fitted
    met = unique(met);                  % some metabolites may have multiple fragments. this step eliminates duplicates
    M0 = or(M0,ismember([mfamodel.mapnode.met]',met));
end
mmm = full(mfamodel.mmm);
mmm = mmm + eye(size(mmm));         %matrix indicating whether metabolite j is produced from metabolite i directly in any reaction
M = backtrack(mmm,M0);              %identify list of metabolites to be traced
mfamodel.mapnode = mfamodel.mapnode(M);
for i = 1:length(mfamodel.mapnode)
    mfamodel.mapnode(i).ammC = [mfamodel.mapnode(i).ammC(M)]';
    mfamodel.mapnode(i).ammC{i} = eye(size(mfamodel.mapnode(i).ammC{i}));
end


natms = sum([mfamodel.mapnode.nC]);
M0 = false(natms,1);
for nx = 1:nxp
    ndata = length([data(nx).msdata]);
    for i = 1:ndata
        met = data(nx).msdata(i).met;
        mpos = find(ismember([mfamodel.mapnode.met]',met));
        ati = sum([mfamodel.mapnode(1:mpos-1).nC])+1;
        atf = sum([mfamodel.mapnode(1:mpos).nC]);
        m0 = M0(ati:atf);
        [~,atms] = emu2info(data(nx).msdata(i).emu);
        m1 = false(size(m0));
        m1(atms) = true;
        M0(ati:atf) = or(m0,m1);
    end
end

M = backtrack(cell2mat([mfamodel.mapnode.ammC]'),M0);   %identify atoms to be traced

mcut = false(size(mfamodel.mapnode));
start = 0;
for i = 1:length(mcut)
    mcx = M(start+1:(start+mfamodel.mapnode(i).nC));
    mfamodel.mapnode(i).emu1 = mfamodel.mapnode(i).emu1(mcx);
    mcut(i) = any(mcx);
    start = start + mfamodel.mapnode(i).nC;
end


mfamodel.mapnode = mfamodel.mapnode(mcut);
%for i = 1:length(mfamodel.mapnode)
%    mfamodel.mapnode(i).ammC = mfamodel.mapnode(i).ammC(mcut);
%end

% identify EMUs involved in EMU balances
% Enumerating all possible EMUs given involved atoms
nemu = zeros(size(mfamodel.mapnode));
nemu = nemu(:);
for i = 1:length(mfamodel.mapnode)
    %disp(i)
    mfamodel.mapnode(i).emu1 = makeallmu(mfamodel.mapnode(i).emu1);
end
emus = [mfamodel.mapnode.emu1];
emuids = {emus.id};
% correction for symmetry:
nosym = true(length(emuids),1);
mets = [mfamodel.mapnode.met]';
for i = 1:length(emuids)
    m = emu2info(emuids{i});
    if isempty(mfamodel.mapnode(ismember(mets,m)).symm)
        nosym(i) = false;
    end
end
Semu = emuids(nosym);
acfor = false(size(Semu));
symgrp = cell(size(Semu));
for i = 1:length(Semu)
    [m,atms] = emu2info(Semu{i});
    p = ismember(mets,m);
    symm = str2num(mfamodel.mapnode(p).symm);
    satm = sort(symm(atms));
    symgrp{i} = info2emu(m,satm);
    if any(ismember(Semu(1:i),symgrp(i)))
        acfor(i) = true;
    end
end
Semu = Semu(~acfor);
symgrp = symgrp(~acfor);
for c1 = 1:length(mfamodel.mapnode)
    for c2 = length(mfamodel.mapnode(c1).emu1):-1:1
        for c3 = 1:length(mfamodel.mapnode(c1).emu1(c2).esrc)
            for c4 = 1:length(mfamodel.mapnode(c1).emu1(c2).esrc(c3).atsrc)
                p1 = ismember(symgrp,mfamodel.mapnode(c1).emu1(c2).esrc(c3).atsrc(c4));
                if any(p1)
                    mfamodel.mapnode(c1).emu1(c2).esrc(c3).atsrc(c4) = Semu(p1);
                end
            end
        end
        p1 = ismember(symgrp,mfamodel.mapnode(c1).emu1(c2).id);
        if any(p1)
            ex = [mfamodel.mapnode(c1).emu1];
            id = {ex.id};
            p2 = ismember(id,Semu(p1));
            mfamodel.mapnode(c1).emu1(p2).esrc = [mfamodel.mapnode(c1).emu1(p2).esrc,mfamodel.mapnode(c1).emu1(c2).esrc];
            mfamodel.mapnode(c1).emu1(c2) = [];
        end
    end
    nemu(c1) = length(mfamodel.mapnode(c1).emu1);
end
emus = [mfamodel.mapnode.emu1];
emuids = {emus.id};

Se = zeros(sum(nemu));
for i = 1:length(emus)
    if emus(i).bal
    esrc = [emus(i).esrc];
    esc = [esrc.atsrc];
    Se(i,:) = double(ismember(emuids,esc));
    end
end


Se = sparse(Se) + speye(sum(nemu));
M0 = false(length(emuids),1);
for i = 1:nxp
    mes = unique([data(i).msdata.emu]);
    M0 = or(M0,ismember(emuids',mes));
end
M = backtrack(Se,M0);   %identifying set of EMUs to be traced
emuids = emuids(M);
emus = emus(M);
Se = Se(M,:);
Se = Se(:,M);


% constructing A and B matrices
esize = [emus.size]';
[~,q] = sort(esize);
sz = sort(unique(esize));
nsize = length(unique(esize));
emod.mdvsim.X = cell(nsize,nxp);
emod.mdvsim.Y = cell(nsize,nxp);
emod.mdvsim.A = cell(nsize,1);
emod.mdvsim.B = cell(nsize,1);
emod.mdvsim.dAdv = cell(nsize,1);
emod.mdvsim.dBdv = cell(nsize,1);
emod.mdvsim.dAdu = cell(nsize,1);
emod.mdvsim.dBdu = cell(nsize,1);
emod.mdvsim.dcxdc = cell(nsize,1);
emod.mdvsim.cY = cell(nsize,1);
emod.mdvsim.sA = cell(nsize,1);
emod.mdvsim.sB = cell(nsize,1);
emod.mdvsim.sX = cell(nsize,nxp);
emod.mdvsim.sY = cell(nsize,nxp);
emod.mdvsim.Xid = cell(nsize,1);
nemus = zeros(nsize,1);
nip = nemus;
emus = emus(q); %size-sorted emus
balemu = emus([emus.bal]);
esize = [emus.size]';
balsize = [balemu.size]';
%allids = cell(0,1);
N = mfamodel.varmgmt.N;
nu = length(N(1,:));
allmets = [mfamodel.node.met];
for i = 1:nsize
    e1 = emus(ismember(esize,sz(i)));
    bals = [e1.bal];
    ebal = e1(bals);
    ey = e1(~bals);
    nx = length(ebal);
    nemus(i) = nx;
    ny = length(ey);
    nv = length(ebal(1).esrc(1).dAdv);
    %for xpt = 1:nxp
    %    emod.mdvsim.X{i,xpt} = zeros(nx,sz(i)+1);
    %end
    %emod.mdvsim.A0 = sparse(zeros(nx,nx));
    s = 0;
    dadv = cell(nx,nx);
    %dadv = dadu;
    %dadu(1:nx,1:nx) = {sparse(zeros(1,nv)*N)};
    dadv(1:nx,1:nx) = {sparse(zeros(1,nv))};
    dbdvf = cell(nx,ny);
    dbdvf(:,:) = {sparse(zeros(1,nv))};
    %dbduf = cell(nx,ny);
    %dbduf(:,:) = {sparse(zeros(1,nv)*N)};
    %dbduc = cell(nx,0);
    dbdvc = cell(nx,0);
    dcdc = zeros(nx,length(allmets));
    idx = {ebal.id}';
    idy = {ey.id}';
    idm = {e1.id}';
    emod.mdvsim.Xid{i} = idx;
    c = zeros(0,length(balsize)-sum(ismember(balsize,[sz(i):max(esize)])));
    %dcdu = cell(0,1);
    dcdv = cell(0,1);
    jinds = zeros(0,1);
    for j = 1:nx
        met = ebal(j).met;
        dcdc(j,ismember(allmets,met)) = 1;
        for k = 1:length(ebal(j).esrc)
            if ebal(j).esrc(k).isconv
                cy = cell(1,i-1);
                for k1 = 1:i-1
                    cy{1,k1} = zeros(1,length(emod.mdvsim.Xid{k1}));
                end
                for k1 = 1:length(ebal(j).esrc(k).atsrc)
                    for k2 = 1:i-1
                        p = ismember(emod.mdvsim.Xid{k2},ebal(j).esrc(k).atsrc(k1));
                        if any(p)
                            cy{1,k2} = cy{1,k2} + double(p');
                        end
                    end
                end
                c = [c;cell2mat(cy)];
                %dcdu = [dcdu;{sparse(ebal(j).esrc(k).dAdv*N)}];
                dcdv = [dcdv;{sparse(ebal(j).esrc(k).dAdv)}];
                jinds = [jinds;j];
            elseif ismember(ebal(j).esrc(k).atsrc,emod.mdvsim.Xid{i})
                p = ismember(emod.mdvsim.Xid{i},ebal(j).esrc(k).atsrc);
                %dadu{j,p} = dadu{j,p} + sparse(ebal(j).esrc(k).dAdv*N);
                dadv{j,p} = dadv{j,p} + sparse(ebal(j).esrc(k).dAdv);
            else
                p = ismember(idy,ebal(j).esrc(k).atsrc);
                dbdvf{j,p} = dbdvf{j,p} + sparse(ebal(j).esrc(k).dAdv);
                %dbduf{j,p} = dbduf{j,p} + sparse(ebal(j).esrc(k).dAdv*N);
            end
        end
    end
    %Handling fixed term inputs
    for xpt = 1:nxp
        emod.mdvsim.Y{i,xpt} = srcmdv(ey,data(xpt));
        emod.mdvsim.sY{i,xpt} = sparse(zeros(length(ey),nu*(sz(i)+1)));
        emod.mdvsim.X{i,xpt} = sparse(zeros(nx,sz(i)+1));
        emod.mdvsim.sX{i,xpt} = sparse(zeros(nx,nu*(sz(i)+1)));
    end
    
    % handling EMU convolutions
    nconv = size(c,1);
    if nconv > 0
        C1 = c;
        for nc = nconv:-1:1
            p = ismember(C1,C1(nc,:),'rows');
            if sum(p)>1
                C1(nc,:) = [];
            end
        end
        nconv = length(C1(:,1));
        dbdvc = cell(nx,nconv);
        dbdvc(:,:) = {sparse(zeros(1,nv))};
        %dbduc = cell(nx,nconv);
        %dbduc(:,:) = {sparse(zeros(1,nv)*N)};
        for nc = 1:length(jinds)
            p = ismember(C1,c(nc,:),'rows');
            %dbduc{jinds(nc),p} = dbduc{jinds(nc),p} + dcdu{nc};
            dbdvc{jinds(nc),p} = dbdvc{jinds(nc),p} + dcdv{nc};
        end
        cY = [zeros(ny,length(C1(1,:)));C1];
        Yconv = cell(nconv,1);
        Yconv(:,:) = {1};
        for xpt = 1:nxp
            emod.mdvsim.Y{i,xpt} = [emod.mdvsim.Y{i,xpt};Yconv];
            emod.mdvsim.sY{i,xpt} = [emod.mdvsim.sY{i,xpt};zeros(nconv,nu*(sz(i)+1))];
        end
        emod.mdvsim.cY{i} = cY;
        %dbdu = [dbduf,dbduc];
        dbdv = [dbdvf,dbdvc];
    else
        %dbdu = dbduf;
        dbdv = dbdvf;
    end
    nyt = ny+nconv;
    for j = 1:nx
        %sa = zeros(1,nv)*N;
        sa = zeros(1,nv);
        for j1 = 1:nx
            %sa = sa+dadu{j,j1};
            sa = sa+dadv{j,j1};
        end
        for j1 = 1:nyt
            %sa = sa+dbdu{j,j1};
            sa = sa+dbdv{j,j1};
        end
        dadv{j,j} = -sa;
    end
    nip(i) = nyt;
    
    emod.mdvsim.dAdv{i} = cell2mat(dadv(:));
    emod.mdvsim.dAdu{i} = emod.mdvsim.dAdv{i}*N;
    emod.mdvsim.A{i} = zeros(nx,nx);
    emod.mdvsim.dBdv{i} = cell2mat(dbdv(:));
    emod.mdvsim.dBdu{i} = emod.mdvsim.dBdv{i}*N;
    emod.mdvsim.B{i} = zeros(nx,nyt);
    dadu = emod.mdvsim.dAdu{i};
    dadu = reshape(dadu,nx,nx*nu);
    dadu = mat2cell(dadu,nx,nx*ones(1,nu));
    dadu = dadu';
    emod.mdvsim.sA{i} = cell2mat(dadu);
    dbdu = emod.mdvsim.dBdu{i};
    dbdu = reshape(dbdu,nx,nyt*nu);
    dbdu = mat2cell(dbdu,nx,nyt*ones(1,nu));
    dbdu = dbdu';
    emod.mdvsim.sB{i} = cell2mat(dbdu);
    
    %pool size associations
    emod.mdvsim.dcxdc{i} = dcdc;
    
end
emod.mdvsim.nemu = nemus;
emod.mdvsim.nip = nip;
for i = 1:nxp
    emod.mdvsim.Y{1,i} = cell2mat(emod.mdvsim.Y{1,i});
end
%storing meta data for organized output file
emod.vardata = mfamodel.varmgmt;
emod.vardata.flxdata = mfamodel.flux;
%emod.vardata.exptdata = mfamodel.data;
emod.vardata.vb = [1e-7*ones(size(emod.vardata.virr)),1e8*ones(size(emod.vardata.virr))];
emod.vardata.nu = length(emod.vardata.N(1,:));

emod.options = mfamodel.options;
%indexing data for collection after forward simulation

[datx(1,1:nxp)] = deal(struct('exptname','e1','flxval',{1},'flxind',1,...
    'msid','1','msval',{[1]},'msind',{[1]},'mswt',{[1]},'mscorr',{[1]}));
emod.data = datx;
%{
emod.data.flxval = cell(length(mfamodel.data.flux),1);
emod.data.flxind = zeros(length(mfamodel.data.flux),1);
emod.data.flxwt = cell(length(mfamodel.data.flux),1);
emod.data.msval = cell(length(mfamodel.data.msdata),1);
emod.data.msind = cell(length(mfamodel.data.msdata),1);
emod.data.mswt = cell(length(mfamodel.data.msdata),1);
emod.data.mscorr = cell(length(mfamodel.data.msdata),1);
%}
nms = zeros(nxp,1);
nh = nxp;
for i = 1:nxp
    %nms(i) = length([mfamodel.data(i).msdata.mdv]);
    msdatasize = size(mfamodel.data(i).msdata);
    nms(i) = msdatasize(2);
    nh(i) = length(mfamodel.data(i).msdata);
    emod.data(i).nh = nh(i);
end
%emod.Eh = zeros(nms,nh);
%emod.data.exptname = mfamodel.data.name;

% indexing flux measurements for collection
for xpt = 1:nxp
    emod.data(xpt).msval = cell(length(mfamodel.data(xpt).msdata),1);
    emod.data(xpt).msind = cell(length(mfamodel.data(xpt).msdata),1);
    emod.data(xpt).mswt = cell(length(mfamodel.data(xpt).msdata),1);
    emod.data(xpt).mscorr = cell(length(mfamodel.data(xpt).msdata),1);
    emod.data(xpt).msid = cell(length(mfamodel.data(xpt).msdata),1);
    emod.data(xpt).flxval = zeros(length(data(xpt).flux),1);
    emod.data(xpt).flxwt = zeros(length(data(xpt).flux),1);
    emod.data(xpt).flxind = zeros(length(data(xpt).flux),1);
    %emod.data(xpt).flxval = [data(xpt).flux.val]';
    %emod.data(xpt).flxwt = [data(xpt).flux.std]';
    flxs = [emod.vardata.flxdata.name]';
    for i = 1:length(data(xpt).flux)
        emod.data(xpt).flxval(i) = data(xpt).flux(i).val;
        emod.data(xpt).flxwt(i) = data(xpt).flux(i).std;
        emod.data(xpt).flxind(i,1) = find(ismember(flxs,data(xpt).flux(i).name));
    end
    emod.data(xpt).exptname = data(xpt).name;
end

%ms data
%lctr = 0;
for xpt = 1:nxp
    for i = 1:length(mfamodel.data(xpt).msdata)
        emod.data(xpt).msid{i,1} = data(xpt).msdata(i).name{1};
        emod.data(xpt).msval{i,1} = mfamodel.data(xpt).msdata(i).mdv;
        emod.data(xpt).mswt{i,1} = mfamodel.data(xpt).msdata(i).std;
        emod.data(xpt).mscorr{i,1} = mfamodel.data(xpt).msdata(i).mdvreg;
        len = length(emod.data(xpt).msval{i});
        [~,atms] = emu2info(mfamodel.data(xpt).msdata(i).emu{1});
        natms = length(atms);
        iind = find(ismember(sz,natms));
        emupos = find(ismember(emod.mdvsim.Xid{iind},mfamodel.data(xpt).msdata(i).emu));
        emod.data(xpt).msind{i,1} = [xpt,iind,emupos,len];
        %emod.Eh(lctr+1:lctr+len,i) = 1;
        %lctr = lctr+len;
    end
end
        
end


function M = backtrack(m,M0)
% helper function for backtracking
done = false;
M = M0;
while ~done
    M1 = any(m(M,:),1);
    if isequal(M1,M)
        M = or(M,M1);
        done = true;
    else
        M = M1;
    end
end
end


%Helper function to extraction name and atom information from emu name
function [name,atom_ids] = emu2info(emu_id)
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


function emu = makeallmu(emu0)
%Helper function to generate a complete list of EMUs given a set of atoms
eb = struct('id',[],'met',[],'nC',[],'size',[],'proxy',[],'dAdv',[],...
    'dBdv',[],'C',[],'esrc',[],'bal',true,'src',false);
l = length(emu0);
atms = zeros(1,l);
for i = 1:l
    [~,atms(i)] = emu2info(emu0(i).id);
end
Ebx = dec2bin(1:power(2,l));
Ebx = Ebx(1:end-1,2:end);
Ebx = fliplr(Ebx);
Ebin = '';
Ebin(1:size(Ebx,1),1:(2*size(Ebx,2)-1)) = ',';
for j = 1:size(Ebx,2)
    Ebin(:,2*j-1) = Ebx(:,j);
end
Ebin = str2num(Ebin);
nemu = size(Ebin,1);
[emu(1:nemu)] = deal(eb);
for i = 1:nemu
    emu(i).bal = emu0(1).bal;
    emu(i).src = emu0(1).src;
    elen = length(find(Ebin(i,:)));
    if elen < 2
        % size 1 EMU
        emu(i).id = emu0(logical(Ebin(i,:))).id;
        emu(i).met = emu0(logical(Ebin(i,:))).met;
        emu(i).nC = emu0(logical(Ebin(i,:))).nC;
        emu(i).size = elen;
        emu(i).esrc = emu0(logical(Ebin(i,:))).esrc;
    else
        eatms = atms(logical(Ebin(i,:)));
        px = ismember(atms,eatms(1));
        id = [emu0(px).id,','];
        for j = 2:length(eatms)
            id = [id,num2str(eatms(j)),','];
        end
        id = id(1:end-1);
        emu(i).id = id;
        emu(i).met = emu0(logical(Ebin(i,:))).met;
        emu(i).nC = emu0(logical(Ebin(i,:))).nC;
        emu(i).size = elen;
        
        emu(i).esrc = emu0(px).esrc;
        nv = length(emu(i).esrc);
        for j = 1:nv
            ea = cell(length(eatms),1);
            mmm = zeros(length(eatms),length(emu(i).esrc(j).rmmm));
            for k = 1:length(eatms)
                ea(k) = emu0(ismember(atms,eatms(k))).esrc(j).atsrc;
                mmm(k,:) = emu0(ismember(atms,eatms(k))).esrc(j).rmmm;
            end
            ey = cell(sum(any(mmm,1)),1);
            mmm = mmm(:,any(mmm,1));
            for ctr = 1:length(ey)
            k = 1;
            ex = ea(logical(mmm(:,ctr)));
            while k<length(ex)
                m1 = emu2info(ex{k});
                [m2,a2] = emu2info(ex{k+1});
                if isequal(m1,m2)
                    %a = [a1,a2];
                    ex{k} = [ex{k},',',num2str(a2)];
                    ex(k+1) = [];
                else
                    k = k+1;
                end
            end
            [m1,a1] = emu2info(ex{1});
            a1 = sort(a1);
            ey{ctr} = [m1,'-'];
            for z = 1:length(a1)
                ey{ctr} = [ey{ctr},num2str(a1(z)),','];
            end
            ey{ctr} = ey{ctr}(1:end-1);
            end
            emu(i).esrc(j).atsrc = ey';
            if length(ey)>1
                emu(i).esrc(j).isconv = true;
            end
        end
    end
end

end


function emuid = info2emu(met,atms)
%helper function to create EMU from metabolite and atom information
emuid = [met,'-'];
atms = sort(atms);
for i = 1:length(atms)
    emuid = [emuid,num2str(atms(i)),','];
end
emuid = emuid(1:end-1);
end


function Y = srcmdv(srcmu,data)
%Helper function to compute MDV for source EMUs
tracer = data.tracer;
Y = cell(length(srcmu),1);
c = [0.989184 0.010816];
trcmets = [tracer.met]';
for i = 1:length(srcmu)
    met = {srcmu(i).met};
    p = ismember(trcmets,met);
    if ~any(p)
        [~,atms] = emu2info(srcmu(i).id);
        y0 = 1;
        for j = 1:length(atms)
            y0 = conv(y0,c);
        end
        Y{i} = y0;
    else
        trc = tracer(p).trc;
        f = [trc.frac];
        y0 = cell(length(f),1);
        y0(:) = {1};
        [~,atms] = emu2info(srcmu(i).id);
        for j = 1:length(y0)
            for k = 1:length(atms)
                y0{j} = conv(y0{j},trc(j).atmdv(atms(k),:));
            end
        end
        y0 = cell2mat(y0);
        f = f(:);
        y0 = f'*y0;
        Y{i} = y0;
    end
end
end