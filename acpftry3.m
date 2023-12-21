%% Set line resistance to zero
MPCr0 = MPC;
MPCr0.bus(:,PD) = (sum(MPC.gen(:,PG))/sum(MPC.bus(:,PD)))*MPC.bus(:,PD);
MPCr0.branch(:,BR_R) = 0;
MPCr0.bus(:,QD) = 0;
%MPCr0.bus(:,VA) = 0; %Can comment out this line because it converges flat start too
MPCr0 = runpf(MPCr0);

%% Increase resistance by step
nsteps = 200;
Rstep = (1/nsteps)*MPC.branch(:,BR_R);
SLACKGEN = MPCr0.gen(:,GEN_BUS)==MPCr0.bus(MPCr0.bus(:,BUS_TYPE)==REF, BUS_I);
VMlog = zeros(nbus, nsteps);
VAlog = zeros(nbus, nsteps);
slacklog = zeros(1, nsteps);
slack1log = zeros(1, nsteps);

Pold = 0;
for i = 1:nsteps
MPCr0.branch(:,BR_R) = MPCr0.branch(:,BR_R) + Rstep;
MPCr1 = runpf(MPCr0);
Pnew = sum(abs(MPCr1.branch(:,PT)+MPCr1.branch(:,PF)));
MPCr0.bus(:,PD) = MPCr0.bus(:,PD) - ...
    ((Pnew-Pold)/sum(MPCr0.bus(:,PD)))*MPCr0.bus(:,PD);

% This converges
%MPCr0.bus(:,PD) = MPCr0.bus(:,PD) - ...
%    0.02*(sum(abs(MPCr1.branch(:,PT)+MPCr1.branch(:,PF)))/sum(MPCr0.bus(:,PD)))*MPCr0.bus(:,PD);

MPCr0 = runpf(MPCr0);
Pold = Pnew;

VMlog(:,i) = MPCr0.bus(:,VM);
VAlog(:,i) = MPCr0.bus(:,VA);
slacklog(i) = sum(MPCr0.gen(SLACKGEN,PG) - MPC.gen(SLACKGEN,PG));
slack1log(i) = sum(MPCr1.gen(SLACKGEN,PG) - MPC.gen(SLACKGEN,PG));

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

for i = 1:5:nconverged
    subplot(1,2,1);
    title(sprintf('Voltage RMS, %d/%d', i, nconverged));
    scatter(node_tbl.LONGITUDE, node_tbl.LATITUDE, 5, VMlog(:,i), 'filled');

    
    subplot(1,2,2);
    title(sprintf('Voltage Angle, %d/%d', i, nconverged));
    scatter(node_tbl.LONGITUDE, node_tbl.LATITUDE, 5, VAlog(:,i), 'filled');

    pause(0.0001);
end
end