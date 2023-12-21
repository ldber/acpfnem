%% Set line resistance to zero
MPCr0 = MPCdc;
MPCr0.branch(:,BR_R) = 0;
MPCr0.bus(:,QD) = 0;
%MPCr0.bus(:,VA) = 0; %Can comment out this line because it converges flat start too
MPCr0 = runpf(MPCr0);

%% Increase resistance by step
Rstep = 0.01*MPCdc.branch(:,BR_R);
SLACKGEN = MPCr0.gen(:,GEN_BUS)==MPCr0.bus(MPCr0.bus(:,BUS_TYPE)==REF, BUS_I);
VMlog = zeros(nbus, 100);
VAlog = zeros(nbus, 100);
slacklog = zeros(1, 100);
slack1log = zeros(1, 100);


for i = 1:100
MPCr0.branch(:,BR_R) = MPCr0.branch(:,BR_R) + Rstep;
MPCr1 = runpf(MPCr0);
MPCr0.bus(:,PD) = MPCr0.bus(:,PD) - ...
    (sum(abs(MPCr1.branch(:,PT)+MPCr1.branch(:,PF)))/sum(MPCr0.bus(:,PD)))*MPCr0.bus(:,PD);
MPCr0 = runpf(MPCr0);

VMlog(:,i) = MPCr0.bus(:,VM);
VAlog(:,i) = MPCr0.bus(:,VA);
slacklog(i) = sum(MPCr0.gen(SLACKGEN,PG) - MPCdc.gen(SLACKGEN,PG));
slack1log(i) = sum(MPCr1.gen(SLACKGEN,PG) - MPCdc.gen(SLACKGEN,PG));

if ~MPCr0.success
    break
end
end

%% Plotting Slack Generation
figure;
plot(1:i, slacklog(1:i), 'b', 'DisplayName', 'After Demand Adjustment'); hold on;
plot(1:i, slack1log(1:i), 'r', 'DisplayName', 'Before Demand Adjustment');