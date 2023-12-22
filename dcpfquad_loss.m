% It has been observed that the power flow converges when resistances of
% the lines are set to zero. However, when the line resistances are brought
% to their nominal values, the nodal demands need to be adjusted for the AC
% power flow to converge. Rather that performing this adjustment in an ad
% hoc, trial and error way, optimisation can be used. This script is a
% preliminary investigation.

%% Initialise to power flow that converges
% The demands in MPCdc were calculated in dcpfquad, and are optimised for
% no line losses. Since this power flow converges, it is used as the
% starting point.

MPCr0 = MPCdc;
MPCr0.branch(:,BR_R) = 0;
MPCr0.bus(:,QD) = 0;
MPCr0 = runpf(MPCr0);

%% Add a bit of line resistance
% If all the resistance is included at once, the power flow does not
% converge (to see this, set nsteps=1 below). Hence the resistances are
% incremented by a small amount. Then a power flow can be run. Since the
% network has not changed much, it should converge. From the results of the
% power flow, the active power lost in each line is found. 

nsteps = 200;
Rstep = (1/nsteps)*MPCdc.branch(:,BR_R);
MPCr0.branch(:,BR_R) = MPCr0.branch(:,BR_R) + Rstep;
MPCr1 = runpf(MPCr0);
Plineloss = abs(MPCr1.branch(:,PT)+MPCr1.branch(:,PF));

% These losses are being supplied by the slack generator which causes its
% dispatch to deviate from the historic dispatch. Thus, the nodal demands
% need to be adjusted for a better solution.


%% Defining Optimisation Objective
% As in dcpfquad, (delta_k-delta_l)/x_kl is approximately the active power
% that flows from bus k to bus l. This is conveniently expressed int matrix
% form using the incidence matrix: Pflow = D*I*Vang.

Incidence = zeros(nline, nbus);
for i = 1:nline
    Incidence(i, node_tbl.NODE_ID == edge_tbl(i,:).FROM_NODE) = -1;
    Incidence(i, node_tbl.NODE_ID == edge_tbl(i,:).TO_NODE) = 1;
end

D = diag(1./((edge_tbl.X_PU)./edge_tbl.NUM_LINES));

% In dcpfquad, it is assumed there are no transmission losses, hence the
% sending power is equal to the receiving power. Thus, the power
% 'accumulating' at each bus could be easily calculated by premultiplying
% the power flows with the transpose of the incidence matrix. However,
% since the system now contains losses, these are accounted for by creating
% a modified incidence matrix.

IncidenceLoss = zeros(nline, nbus);
for i = 1:nline
    if MPCr1.branch(i,PF) > 0
        IncidenceLoss(i, MPCr1.bus(:,BUS_I)==MPCr1.branch(i,F_BUS)) = -1;
        IncidenceLoss(i, MPCr1.bus(:,BUS_I)==MPCr1.branch(i,T_BUS)) = -MPCr1.branch(i,PT)/MPCr1.branch(i,PF);
    else
        IncidenceLoss(i, MPCr1.bus(:,BUS_I)==MPCr1.branch(i,F_BUS)) = MPCr1.branch(i,PF)/MPCr1.branch(i,PT);
        IncidenceLoss(i, MPCr1.bus(:,BUS_I)==MPCr1.branch(i,T_BUS)) = 1;        
    end
end

% Now, the power 'accumulating' at each bus can be written as:
% Pbus=IncidenceLoss'*D*Incidence*Vang. To minimise the angle differences
% in the system, the objective function is set up as follows.

IiIltDI = Incidence*pinv(IncidenceLoss'*D*Incidence);
Hl = 2*(IiIltDI')*IiIltDI;
f = zeros(nbus, 1);

%% Defining Constraint
% The same constraints are used as in dcpfquad. Namely, the generator
% dispatches are +/- 1MW from their historic dispatch and the demands are
% constrained to +/- 10%.

A = [ones(1,nbus); -ones(1,nbus)];
b = [1e-04; 1e-04];

Aeq = [];
beq = [];

LB = zeros(nbus,1);
UB = zeros(nbus,1);

for i = 1:(ngen+nhvdc)
    mask = MPC.bus(:,BUS_I) == MPC.gen(i,GEN_BUS);
    LB(mask) = LB(mask) + MPC.gen(i,PG) - 1;
    UB(mask) = UB(mask) + MPC.gen(i,PG) + 1;
end

for i = 1:nbus
    if MPCdc.bus(i,PD) > 0
        LB(i) = LB(i) - 1.1*MPCr1.bus(i,PD);
        UB(i) = UB(i) - 0.9*MPCr1.bus(i,PD);
    else
        LB(i) = LB(i) - 0.9*MPCr1.bus(i,PD);
        UB(i) = UB(i) - 1.1*MPCr1.bus(i,PD);        
    end
end

%% Running Optimisation
Popt = quadprog(Hl,f,A,b,Aeq,beq,LB,UB);

%% Updating Demands
for i = 1:nbus
    mask = MPC.gen(:,GEN_BUS) == MPC.bus(i,BUS_I);
    MPCr0.bus(i,PD) = median([MPCr1.bus(i,PD)*[0.9,1.1], sum(MPC.gen(mask,PG)) - Popt(i)]);
    if any(mask)
        MPCr0.gen(mask,PG) = MPC.gen(mask,PG) + ...
                (Popt(i)-sum(MPC.gen(mask,PG))+MPCr0.bus(i,PD))/sum(mask);
    end
end

%% Running Power Flow
MPCr0 = runpf(MPCr0);
