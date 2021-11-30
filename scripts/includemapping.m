function [ model,mfamodel ] = includemapping( model,filename )
%INCLUDEMAPPING upgrades the existing metabolic model to include atom
%transition data the upgraded model separates reversible reactions into
%forward and backward reactions.In addition, the updated model also lists
%the sources for all mapped metabolites in order to simplify path tracing
%for MFA. This function upgrades the model from level 1 to level 2.
%   [MODEL,MFAMODEL] = includemapping(MODEL,FILENAME) reads reaction atom 
%   mapping information contained within the excel file FILENAME and 
%   incorporates into the RXNPROP field of MODEL. A subfield (FLUX) is
%   added to RXNPROP summarizing key information necessary for carbon path 
%   tracing. A new substructure called NODES is added to MODEL describing 
%   carbon atom sources for all nodes contained in METPROP to simplify EMU 
%   network generation.
%   MFAMODEL is a modified network structure containing properties
%   pertinent to 13C-MFA. It contains the following fields:
%       FLUX:   All reactions decomposed into forward and reversible
%               reactions with *.f representing th forward reaction and *.b
%               representing the reverse reaction. This process is to
%               ensure non-negativity of all estimated fluxes.
%       NODE:   Contains the name of every metabolite, with indicators of
%               whether it is balanced, internal/source/sink metabolite,
%               and whether it is involved involved in the overall atom
%               mapping network.
%       MAPNODE:Contains details about symmetry and the metabolites that
%               are the direct precursors of this particular metabolite.
%               This data is used to create a genome-scale connectivity
%               matrix. The precursor to the size-1 EMU network is also
%               generated.
%       MMM:    This is the genome-scale connectivity matrix with row
%               i containing the list of precursor metabolites, and colunm
%               j containing the list of products produced from metabolite
%               j. It is a rectangular matrix with source metabolites
%               excluded from the set of metabolites forming the rows.
%       DATA:   Will be described in Level 3
%       STRC:   Structural parameters including the EMU network, described
%               Level 4.
%       LEVEL:  Current Level is 2, indicating that it can be used for
%               metabolite path tracing.
%
%   Written by Saratram Gopalakrishnan on 09/03/2016 at 102 MRL, Penn State

%initialization
mapstr = struct('reactants',[],'products',[],'amm',[],'mmm',[]);
flx = struct('type',[],'name',[],'reactant',[],'rstoic',[],'product',[],'pstoic',[],'mapping',mapstr,'main',[],'dilution',[]);
nodestr = struct('met',[],'bal',[],'type',[],'map',[]);
nmapstr = struct('met',[],'defm',[],'symm',[],'dmdv',[],'bal',[],'src',[],'nC',[],'ammC',[]);
mfamodel = struct('flux',[],'node',nodestr,'mapnode',nmapstr,'data',[],'strc',[],'S',[],'level',2);
source = struct('flx',[],'atsrc',[],'isconv',false,'rmmm',[],'dAdv',[]);
emu = struct('id',[],'met',[],'nC',[],'size',[],'proxy',[],'dAdv',[],'dBdv',[],'C',[],'esrc',[],'bal',true,'src',false);

S_full = zeros(size(model.S,1),size(model.S,2)+sum([model.rxnprop.rev]));
ctr = 0;
virr = false(size(model.S,2)+sum([model.rxnprop.rev]),1);
vfwd = virr;
vrev = virr;
vnetirr = false(size(model.S,2),1);
vnetrev = vnetirr;
%splitting reactions into forward and reverse
for i = 1:length(model.rid)
    if model.rxnprop(i).rev
        vnetrev(i) = true;
        f1 = flx;
        f1.type = 'fwd';
        f1.name = {[char(model.rid(i)),'.f']};
        f1.reactant = model.rxnprop(i).baseflx.reactant;
        f1.rstoic = model.rxnprop(i).baseflx.rstoic;
        f1.product = model.rxnprop(i).baseflx.product;
        f1.pstoic = model.rxnprop(i).baseflx.pstoic;
        f1.main = ~ismember(model.rxnprop(i).subsystem,{'dilution'});
        f1.dilution = ismember(model.rxnprop(i).subsystem,{'dilution'});
        f2 = flx;
        ctr = ctr+1;
        S_full(:,ctr) = model.S(:,i);
        vfwd(ctr) = true;
        f2.type = 'rev';
        f2.name = {[char(model.rid(i)),'.r']};
        f2.reactant = model.rxnprop(i).baseflx.product;
        f2.rstoic = model.rxnprop(i).baseflx.pstoic;
        f2.product = model.rxnprop(i).baseflx.reactant;
        f2.pstoic = model.rxnprop(i).baseflx.rstoic;
        f2.main = f1.main;
        f2.dilution = f1.dilution;
        model.rxnprop(i).flux = [f1,f2];
        ctr = ctr+1;
        S_full(:,ctr) = -S_full(:,ctr-1);
        vrev(ctr) = true;
    else
        f1 = flx;
        f1.type = 'irr';
        f1.name = {[char(model.rid(i)),'.f']};
        f1.reactant = model.rxnprop(i).baseflx.reactant;
        f1.rstoic = model.rxnprop(i).baseflx.rstoic;
        f1.product = model.rxnprop(i).baseflx.product;
        f1.pstoic = model.rxnprop(i).baseflx.pstoic;
        f1.main = ~ismember(model.rxnprop(i).subsystem,{'dilution'});
        f1.dilution = ismember(model.rxnprop(i).subsystem,{'dilution'});
        model.rxnprop(i).flux = f1;
        ctr = ctr+1;
        S_full(:,ctr) = model.S(:,i);
        virr(ctr) = true;
        vnetirr(i) = true;
    end
end
ubf = [model.metprop.ubflag]';
src = [model.metprop.src]';
snk = [model.metprop.snk]';
bal = true(size(ubf));
bal(ubf) = false;
bal(src) = false;
bal(snk) = false;
S_bal = S_full(bal,:);



%extracting raw data
[~,~,raw1] = xlsread(filename,'Atom Transition');
[~,~,raw2] = xlsread(filename,'Symmetry');

%updating symmetry 
mets = raw2(2:end,1);
sym = raw2(2:end,3);
for i = 1:length(mets)
    metlist = [model.metprop.metid]';
    pos = ismember(metlist,mets(i));
    model.metprop(pos,1).symmap = sym{i};
end

%including atom transition information
ctrns = ismember(raw1(:,2),'C');
atmps = raw1(ctrns,:);
[model.metprop(1:end).mapping] = deal(false);
mts = unique(atmps(:,3));
[model.metprop(ismember(metlist,mts)).mapping] = deal(true);
rxns = unique(atmps(:,6));
relmap = ismember(rxns,model.rid);
rxns = rxns(relmap);
atmps = atmps(ismember(atmps(:,6),rxns),:);
chkmap = false(size(model.rid));
for i = 1:1:length(rxns)
    %i;
    %if i == 56
    %    i
    %end
    ptn = ismember(model.rid,rxns(i));
    tmp = atmps(ismember(atmps(:,6),rxns(i)),:);
    m1 = mapstr;
    % reactants
    t1 = tmp(ismember(tmp(:,4),'reactant'),:);
    m1.reactants.name = t1(:,3);
    for j = 1:length(t1(:,1))
        if isnumeric(t1{j,1})
            t1{j,1} = num2str(t1{j,1});
        end
        m1.reactants.tr{j} = strsplit(char(t1(j,1)),',');
    end
    %products
    t1 = tmp(ismember(tmp(:,4),'product'),:);
    m1.products.name = t1(:,3);
    for j = 1:length(t1(:,1))
        if isnumeric(t1{j,1})
            t1{j,1} = num2str(t1{j,1});
        end
        m1.products.tr{j} = strsplit(char(t1(j,1)),',');
    end
    %amm
    m1.amm = cell(length(m1.products),length(m1.reactants));
    for c1 = 1:length(m1.products.name)
        l1 = length(m1.products.tr{c1});
        z1 = ismember([model.metprop.metid],m1.products.name(c1));
        model.metprop(z1).nC = max(model.metprop(z1).nC,l1);
      for c2 = 1:length(m1.reactants.name)
          l2 = length(m1.reactants.tr{c2});
          z1 = ismember([model.metprop.metid],m1.reactants.name(c2));
          model.metprop(z1).nC = max(model.metprop(z1).nC,l2);
          m1.amm{c1,c2} = sparse(zeros(l1,l2));
          for k = 1:l1
              ismp = ismember(m1.products.tr{c1}(k),m1.reactants.tr{c2});
              if ismp
                  [ndx] = ismember(m1.reactants.tr{c2},m1.products.tr{c1}(k));
                  m1.amm{c1,c2}(k,ndx) = 1;
                  m1.mmm(c1,c2) = 1;
              end
          end
      end
    end
          m2 = m1;
          m2.reactants = m1.products;
          m2.products = m1.reactants;
          m2.amm = m2.amm';
          m2.mmm = m2.mmm';
          [l,r] = size(m2.amm);
          for c1 = 1:l
              for c2 = 1:r
                  m2.amm{c1,c2} = m2.amm{c1,c2}';
              end
          end
          if sum(ismember(m1.reactants.name,model.rxnprop(ptn).flux(1).reactant)) == length(m1.reactants.name)
              if model.rxnprop(ptn).rev
                  model.rxnprop(ptn).flux(1).mapping = m1;
                  model.rxnprop(ptn).flux(2).mapping = m2;
              else
                  model.rxnprop(ptn).flux.mapping = m1;
              end
          elseif sum(ismember(m1.reactants.name,model.rxnprop(ptn).flux(1).product)) == length(m1.reactants.name)
              if model.rxnprop(ptn).rev
                  model.rxnprop(ptn).flux(1).mapping = m2;
                  model.rxnprop(ptn).flux(2).mapping = m1;
              else
                  model.rxnprop(ptn).flux.mapping = m2;
              end
          else
              chkmap(ptn) = true;
          end
end
      
if any(chkmap)
disp(['check mapping for ',num2str(find(full(chkmap')))])
end
model.level = 2;
%genome-scale metabolite mapping matrix

%reactions
%null space matrix:
[R,jb] = rref(S_bal);
vfree = true(size(virr));
vfree(jb) = false;
N1 = S_bal(:,vfree);
N2 = S_bal(:,~vfree);
Ndep = -N2\N1;
N = zeros(length(vfree),sum(vfree));
%N(~vfree,:) = -R(1:length(jb),vfree);
N(~vfree,:) = Ndep;
N(vfree,:) = eye(sum(vfree));


mfamodel.flux = [model.rxnprop.flux];
mfamodel.varmgmt.S = S_full;
mfamodel.varmgmt.S_bal = S_bal;
mfamodel.varmgmt.vnetirr = vnetirr;
mfamodel.varmgmt.vnetrev = vnetrev;
mfamodel.varmgmt.virr = virr;
mfamodel.varmgmt.vfwd = vfwd;
mfamodel.varmgmt.vrev = vrev;
mfamodel.varmgmt.N = N;
mfamodel.varmgmt.vfree = vfree;
fluxes = mfamodel.flux;
%metabolites
[mfamodel.node(1:length(model.metprop))] = deal(nodestr);
k = [model.metprop.ubflag]' | [model.metprop.src]' | [model.metprop.snk]';
for i = 1:length(model.metprop)
    mfamodel.node(i).met = model.metprop(i).metid;
    mfamodel.node(i).bal = ~k(i);
    if model.metprop(i).src
        mfamodel.node(i).type = 'src';
    elseif model.metprop(i).snk
        mfamodel.node(i).type = 'snk';
    else
        mfamodel.node(i).type = 'int';
    end
    mfamodel.node(i).map = model.metprop(i).mapping;
end

%mapped metabolites
k = [mfamodel.node.map]';
[mfamodel.mapnode(1:sum(k))] = deal(nmapstr);
mn = mfamodel.node(k);
for i = 1:sum(k)
    mfamodel.mapnode(i).met = mn(i).met;
    mfamodel.mapnode(i).bal = mn(i).bal;
    mfamodel.mapnode(i).src = ismember(mn(i).type,{'src'});
    mfamodel.mapnode(i).symm = model.metprop(ismember([model.metprop.metid],mfamodel.mapnode(i).met)).symmap;
    mfamodel.mapnode(i).dmdv = sparse(zeros(length(mfamodel.flux),sum(k)));
    mfamodel.mapnode(i).nC = model.metprop(ismember([model.metprop.metid],mfamodel.mapnode(i).met)).nC;
    mfamodel.mapnode(i).ammC = cell(1,sum(k));
    ex = emu;
    [ex(1:mfamodel.mapnode(i).nC)] = deal(emu);
    mfamodel.mapnode(i).emu1 = ex;
    for j = 1:mfamodel.mapnode(i).nC
        mfamodel.mapnode(i).emu1(j).met = char(mfamodel.mapnode(i).met);
        mfamodel.mapnode(i).emu1(j).id = [mfamodel.mapnode(i).emu1(j).met,'-',num2str(j)];
        mfamodel.mapnode(i).emu1(j).size = 1;
        mfamodel.mapnode(i).emu1(j).nC = mfamodel.mapnode(i).nC;
        mfamodel.mapnode(i).emu1(j).bal = mfamodel.mapnode(i).bal;
        mfamodel.mapnode(i).emu1(j).src = mfamodel.mapnode(i).src;
    end
end
srt = [mfamodel.mapnode.src];
[~,q] = sort(srt);
mfamodel.mapnode = mfamodel.mapnode(q);
for i = 1:length(mfamodel.mapnode)
    %if ~mfamodel.mapnode(i).src
    for j = 1:length(mfamodel.mapnode)
        mfamodel.mapnode(i).ammC{1,j} = sparse(zeros(mfamodel.mapnode(i).nC,mfamodel.mapnode(j).nC));
    end
    %end
end

for i = 1:length(fluxes)
    if ~isempty(fluxes(i).mapping.reactants)
        for j1 = 1:length(fluxes(i).mapping.products.name)
            for j2 = 1:length(fluxes(i).mapping.reactants.name)
                if any(any(fluxes(i).mapping.amm{j1,j2}))
                    k1 = ismember([mfamodel.mapnode.met],fluxes(i).mapping.products.name(j1));
                    k2 = ismember([mfamodel.mapnode.met],fluxes(i).mapping.reactants.name(j2));
                    mfamodel.mapnode(k1).dmdv(i,k2) = 1;
                    mfamodel.mapnode(k1).ammC{1,k2} = double(or(logical(mfamodel.mapnode(k1).ammC{1,k2}),fluxes(i).mapping.amm{j1,j2}));
                    tamm = fluxes(i).mapping.amm{j1,j2};
                    z1 = any(tamm,2);
                    z2 = 1:length(tamm(1,:));
                    for w = 1:length(z1)
                        if z1(w)
                           e = source;
                           e.flx = fluxes(i).name;
                           e.atsrc = {[char(fluxes(i).mapping.reactants.name{j2}),'-',num2str(z2(logical(tamm(w,:))))]};
                           e.rmmm = zeros(1,length(fluxes(i).mapping.reactants.name));
                           e.rmmm(j2) = 1;
                           e.dAdv = sparse(zeros(1,length(fluxes)));
                           e.dAdv(i) = 1;
                           mfamodel.mapnode(k1).emu1(w).esrc = [mfamodel.mapnode(k1).emu1(w).esrc,e];
                        end
                    end

                end
            end
        end
    end
end

%Metabolite mapping matrix
gsmmm = reshape(any([mfamodel.mapnode.dmdv],1),[sum(k),sum(k)]);
mfamodel.mmm = gsmmm';
end


