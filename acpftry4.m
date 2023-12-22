%% Set line resistance to zero
MPCr0 = MPCdc;
MPCr0.branch(:,BR_R) = 0;
MPCr0.bus(:,QD) = 0;
MPCr0 = runpf(MPCr0);

%% Define Incidence Matrix
Incidence = zeros(nline, nbus);
for i = 1:nline
    Incidence(i, node_tbl.NODE_ID == edge_tbl(i,:).FROM_NODE) = -1;
    Incidence(i, node_tbl.NODE_ID == edge_tbl(i,:).TO_NODE) = 1;
end

D = diag(1./((edge_tbl.X_PU)./edge_tbl.NUM_LINES));

A = [ones(1,nbus); -ones(1,nbus)];
b = [1e-04; 1e-04];

Aeq = [];
beq = [];

%% Increase resistance by step
nsteps = 100;
Rstep = (1/nsteps)*MPCdc.branch(:,BR_R);
SLACKGEN = MPCr0.gen(:,GEN_BUS)==MPCr0.bus(MPCr0.bus(:,BUS_TYPE)==REF, BUS_I);
VMlog = zeros(nbus, nsteps);
VAlog = zeros(nbus, nsteps);
slacklog = zeros(1, nsteps);
slack1log = zeros(1, nsteps);

for i = 1:nsteps
    MPCr0.branch(:,BR_R) = MPCr0.branch(:,BR_R) + Rstep;
    MPCr1 = runpf(MPCr0);

    IncidenceLoss = zeros(nline, nbus);
    for k = 1:nline
        if MPCr1.branch(k,PF) > 0
            IncidenceLoss(k, MPCr1.bus(:,BUS_I)==MPCr1.branch(k,F_BUS)) = -1;
            IncidenceLoss(k, MPCr1.bus(:,BUS_I)==MPCr1.branch(k,T_BUS)) = -MPCr1.branch(k,PT)/MPCr1.branch(k,PF);
        else
            IncidenceLoss(k, MPCr1.bus(:,BUS_I)==MPCr1.branch(k,F_BUS)) = MPCr1.branch(k,PF)/MPCr1.branch(k,PT);
            IncidenceLoss(k, MPCr1.bus(:,BUS_I)==MPCr1.branch(k,T_BUS)) = 1;        
        end
    end

    IiIltDI = Incidence*pinv(IncidenceLoss'*D*Incidence);
    Hl = 2*(IiIltDI')*IiIltDI;
    f = zeros(nbus, 1);

    LB = zeros(nbus,1);
    UB = zeros(nbus,1);

    for k = 1:(ngen+nhvdc)
        mask = MPC.bus(:,BUS_I) == MPC.gen(k,GEN_BUS);
        LB(mask) = LB(mask) + MPC.gen(k,PG) - 1;
        UB(mask) = UB(mask) + MPC.gen(k,PG) + 1;
    end

    for k = 1:nbus
        if MPCdc.bus(k,PD) > 0
            LB(k) = LB(k) - 1.1*MPCr1.bus(k,PD);
            UB(k) = UB(k) - 0.9*MPCr1.bus(k,PD);
        else
            LB(k) = LB(k) - 0.9*MPCr1.bus(k,PD);
            UB(k) = UB(k) - 1.1*MPCr1.bus(k,PD);        
        end
    end

    Popt = quadprog(Hl,f,A,b,Aeq,beq,LB,UB);

    for k = 1:nbus
        mask = MPC.gen(:,GEN_BUS) == MPC.bus(k,BUS_I);
        MPCr0.bus(k,PD) = median([MPCr1.bus(k,PD)*[0.9,1.1], sum(MPC.gen(mask,PG)) - Popt(k)]);
        if any(mask)
            MPCr0.gen(mask,PG) = MPC.gen(mask,PG) + ...
                    (Popt(k)-sum(MPC.gen(mask,PG))+MPCr0.bus(k,PD))/sum(mask);
        end
    end


    MPCr0 = runpf(MPCr0);

    VMlog(:,i) = MPCr0.bus(:,VM);
    VAlog(:,i) = MPCr0.bus(:,VA);
    slacklog(i) = sum(MPCr0.gen(SLACKGEN,PG) - MPCdc.gen(SLACKGEN,PG));
    slack1log(i) = sum(MPCr1.gen(SLACKGEN,PG) - MPCdc.gen(SLACKGEN,PG));

    if ~MPCr0.success
        break
    end
end

nconverged = i;

%% Plotting Slack Generation
figure;
plot((1:nconverged)/nsteps*100, slacklog(1:nconverged), 'b', 'DisplayName', 'After Demand Adjustment'); hold on;
plot((1:nconverged)/nsteps*100, slack1log(1:nconverged), 'r', 'DisplayName', 'Before Demand Adjustment');
xlim([0,100]);
legend();

%% Plotting Voltage RMS and Angle
if 1
aus = shaperead("aus.shp");
fig = figure;
fig.WindowState = 'maximized';

for i = 1:2
    subplot(1,2,i);
    plot(aus(1).X, aus(1).Y, 'k'); hold on;
    plot(aus(2).X, aus(2).Y, 'k');
    plot(aus(3).X, aus(3).Y, 'k');
    plot(aus(4).X, aus(4).Y, 'k');
    axis('square');
    for k = 1:nline
        plot(edge_geotbl.Longitude{k}, edge_geotbl.Latitude{k}, 'Color', [0.7 0.7 0.7]);
    end
end

subplot(1,2,1);
colormap('jet');
caxis([min(min(VMlog(:,1:nconverged))), max(max(VMlog(:,1:nconverged)))]);
colorbar;

subplot(1,2,2);
colormap('jet');
caxis([min(min(VAlog(:,1:nconverged))), max(max(VAlog(:,1:nconverged)))]);
colorbar;

for i = 1:1:nconverged
    subplot(1,2,1);
    title(sprintf('Voltage RMS, %d/%d', i, nconverged));
    scatter(node_tbl.LONGITUDE, node_tbl.LATITUDE, 5, VMlog(:,i), 'filled');

    
    subplot(1,2,2);
    title(sprintf('Voltage Angle, %d/%d', i, nconverged));
    scatter(node_tbl.LONGITUDE, node_tbl.LATITUDE, 5, VAlog(:,i), 'filled');

    pause(0.0001);
end
end