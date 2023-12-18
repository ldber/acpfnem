clc; clear all; close all;

%% Read csv files

node_tbl = readtable('network_nodes.csv');
edge_tbl = readtable('network_edges.csv');
hvdc_tbl = readtable('network_hvdc_links.csv');
gen_tbl = readtable('generators.csv');

nline = size(edge_tbl,1);

nodesTAS = strcmp(node_tbl.STATE_NAME, 'Tasmania');
gensTAS = strcmp(gen_tbl.NEM_REGION, 'TAS1');
linesTAS = false(nline, 1);

for i = 1:nline
    linesTAS(i) = nodesTAS(node_tbl.NODE_ID == edge_tbl(i,:).FROM_NODE) | ...
        nodesTAS(node_tbl.NODE_ID == edge_tbl(i,:).TO_NODE);
end

node_tbl = node_tbl(~nodesTAS, :);
edge_tbl = edge_tbl(~linesTAS, :);
gen_tbl = gen_tbl(~gensTAS, :);

nbus = size(node_tbl,1);
nline = size(edge_tbl,1);
ngen = size(gen_tbl,1);
nhvdc = size(hvdc_tbl,1);

hvdcmainland = false(nhvdc,2);
for i = 1:nhvdc
    hvdcmainland(i,1) = any(hvdc_tbl.FROM_NODE(i)==node_tbl.NODE_ID);
    hvdcmainland(i,2) = any(hvdc_tbl.TO_NODE(i)==node_tbl.NODE_ID);
end
nhvdcbus = sum(sum(hvdcmainland));
hvdc_tbl(strcmp(hvdc_tbl.HVDC_LINK_ID,'DIRECTLINK'),:).HVDC_LINK_ID = {'N-Q-MNSP1'};

NEMREGIONS = {'NSW1';'VIC1';'QLD1';'SA1'};
nreg = size(NEMREGIONS, 1);

%% Make a MATPOWER structure
% To make the MATPOWER structure, we need to specify the base power, bus
% data, branch data and generator data in the required format.

% From Table 2, the base power for per unit is 100MVA
basecase.baseMVA = 100;

% The bus data contains 13 columns. The required format is found using the
% idx_bus function.
basecase.bus = zeros(nbus, 13);
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
basecase.bus(:,BUS_TYPE) = PQ;
for i = 1:ngen
    basecase.bus(node_tbl.NODE_ID == gen_tbl.NODE(i), BUS_TYPE) = PV;
end
basecase.bus(node_tbl.NODE_ID == gen_tbl.NODE(strcmp(gen_tbl.STATIONID, 'MURRAY')), BUS_TYPE) = REF;
basecase.bus(:,BUS_I) = node_tbl.NODE_ID;
basecase.bus(:,BUS_AREA) = int8(categorical(node_tbl.NEM_REGION, NEMREGIONS, 'Ordinal', true));
basecase.bus(:,PD) = node_tbl.PROP_REG_D;
basecase.bus(:,QD) = 0;
basecase.bus(:,VM) = 1;
basecase.bus(:,VA) = 0;
basecase.bus(:,ZONE) = 1;
basecase.bus(:,VMAX) = 1.2;
basecase.bus(:,VMAX) = 0.8;

% The branch data contains 13 columns. The required format is found using
% the idx_brch function.
basecase.branch = zeros(nline, 13);
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
basecase.branch(:,F_BUS) = edge_tbl.FROM_NODE;
basecase.branch(:,T_BUS) = edge_tbl.TO_NODE;
basecase.branch(:,BR_R) = (edge_tbl.R_PU) ./ (edge_tbl.NUM_LINES);
basecase.branch(:,BR_X) = (edge_tbl.X_PU) ./ (edge_tbl.NUM_LINES);
basecase.branch(:,BR_B) = (edge_tbl.B_PU) .* (edge_tbl.NUM_LINES);
basecase.branch(:,BR_STATUS) = 1;
basecase.branch(:,ANGMIN) = -360;
basecase.branch(:,ANGMAX) = 360;

% The generator data contains 21 columns. The required format is found
% using the idx_gen function.
basecase.gen = zeros(ngen+nhvdcbus, 21);
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
basecase.gen(1:ngen, GEN_BUS) = gen_tbl.NODE;
basecase.gen(1:ngen, PMAX) = gen_tbl.REG_CAP;
basecase.gen(1:ngen, PMIN) = -gen_tbl.REG_CAP;
basecase.gen(1:ngen, QMAX) = gen_tbl.REG_CAP;
basecase.gen(1:ngen, QMIN) = -gen_tbl.REG_CAP;
basecase.gen(1:ngen, VG) = 1.02;
basecase.gen(1:ngen, GEN_STATUS) = 1;


% Add generators at start and end nodes of HVDC links
mask = ngen;
for i = 1:nhvdc
    if hvdcmainland(i,1)
        mask = mask+1;
        basecase.bus(hvdc_tbl.FROM_NODE(i)==node_tbl.NODE_ID, BUS_TYPE) = PV;    
        basecase.gen(mask, GEN_BUS) = hvdc_tbl.FROM_NODE(i);
        basecase.gen(mask, PMAX) = hvdc_tbl.REVERSE_LIMIT_MW(i);
        basecase.gen(mask, PMIN) = -hvdc_tbl.FORWARD_LIMIT_MW(i);
        basecase.gen(mask, QMAX) = hvdc_tbl.REVERSE_LIMIT_MW(i);
        basecase.gen(mask, QMIN) = -hvdc_tbl.FORWARD_LIMIT_MW(i);
        basecase.gen(mask, VG) = 1.02;
        basecase.gen(mask, GEN_STATUS) = 1;
    end
    if hvdcmainland(i,2)
        mask = mask+1;
        basecase.bus(hvdc_tbl.TO_NODE(i)==node_tbl.NODE_ID, BUS_TYPE) = PV;    
        basecase.gen(mask, GEN_BUS) = hvdc_tbl.TO_NODE(i);
        basecase.gen(mask, PMAX) = hvdc_tbl.FORWARD_LIMIT_MW(i);
        basecase.gen(mask, PMIN) = -hvdc_tbl.REVERSE_LIMIT_MW(i);
        basecase.gen(mask, QMAX) = hvdc_tbl.FORWARD_LIMIT_MW(i);
        basecase.gen(mask, QMIN) = -hvdc_tbl.REVERSE_LIMIT_MW(i);
        basecase.gen(mask, VG) = 1.02;
        basecase.gen(mask, GEN_STATUS) = 1;
    end
end

%% Read historic demand and dispatch
demand_tbl = readtable('signals_regional_load.csv');
dispatch_tbl = readtable('signals_dispatch.csv', 'PreserveVariableNames', 1);
interconnector_tbl = readtable('interconnector_flows.csv');

timeidx = 24; %This corresponds to 12 noon, 1 June 2017

MPC = loadcase(basecase);

for i = 1:nreg
    mask = MPC.bus(:,BUS_AREA) == i;
    MPC.bus(mask,PD) = MPC.bus(mask,PD)*demand_tbl.(NEMREGIONS{i})(timeidx);
    MPC.bus(mask,QD) = MPC.bus(mask,QD)*demand_tbl.(NEMREGIONS{i})(timeidx);
end

for i = 1:ngen
    mask = strcmp(dispatch_tbl.Properties.VariableNames, gen_tbl.DUID{i});
    if any(mask)
        MPC.gen(i,PG) = dispatch_tbl{timeidx, mask};
    end
end

interconnector_tbl.SETTLEMENTDATE = datetime(interconnector_tbl.SETTLEMENTDATE, 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
timemask = interconnector_tbl.SETTLEMENTDATE == dispatch_tbl.SETTLEMENTDATE(timeidx);

mask = ngen;
for i = 1:nhvdc
    if hvdcmainland(i,1)
        mask = mask + 1;
        MPC.gen(mask,PG) = -interconnector_tbl(timemask&strcmp(interconnector_tbl.INTERCONNECTORID, hvdc_tbl.HVDC_LINK_ID(i)),:).MWFLOW;
    end
    if hvdcmainland(i,2)
        mask = mask + 1;
        MPC.gen(mask,PG) = interconnector_tbl(timemask&strcmp(interconnector_tbl.INTERCONNECTORID, hvdc_tbl.HVDC_LINK_ID(i)),:).MWFLOW;
    end   
end


