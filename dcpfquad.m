% The variables are P at each bus. PQ buses are constrained to their
% historic dispatch. PV buses are constrained to +/- 10% of their historic
% dispatch. The objective is to minimise power flows so that that angle
% differences are minimum hence the DC approximation is more valid. HVDCs
% are set as PQ buses with P corresponding to historical flow data.

% From docs, X = quadprog(H,f,A,b,Aeq,beq,LB,UB)
% Let us specify the H matrix. First find the incidence matrix

Incidence = zeros(nline, nbus);
for i = 1:nline
    Incidence(i, node_tbl.NODE_ID == edge_tbl(i,:).FROM_NODE) = -1;
    Incidence(i, node_tbl.NODE_ID == edge_tbl(i,:).TO_NODE) = 1;
end

D = diag(1./((edge_tbl.X_PU)./edge_tbl.NUM_LINES));

% P = I'*D*I*Vang. Also I*Vang is angle differences. Note that (I'*D) is
% not full rank. The null space of this matrix is spanned by a vector of
% ones. Thus, the angle difference are uniquely specified by
% I*Vang = pinv(I'*D)*P. Now ||I*Vang|| = (I*Vang)'*(I*Vang) =
% (pinv(I'*D)*P)'(pinv(I'*D)*P) = P'*(pinv(I'*D))'*pinv(I'*D)*P. Thus, the
% H matrix for quadprog can be defined as:

iItD = Incidence*pinv(Incidence'*D*Incidence);
H = 2*(iItD')*iItD;

% Note the factor of 2 is because quadprog minimises 0.5*x'*H*x + f'*x. The
% linear term can be set to zero. If considering bidding of generators,
% they can be incorporated into this linear term.

f = zeros(nbus, 1);

% Now to set the constraints. Equality constraints are avoided since they
% might make the optimistation harder to solve. Since the resistors are
% being ignored there are no transmission loss of active power. Thus, it is
% required that sum(P)=0, i.e. the total generation should equal the total
% load. This can be expressed as an inequality constraint:

MPCdc = MPC;
MPCdc.bus(:,PD) = (sum(MPC.gen(:,PG))/sum(MPC.bus(:,PD)))*MPC.bus(:,PD);


A = [ones(1,nbus); -ones(1,nbus)];
b = [1e-04; 1e-04];

% As explained, equality constraints will not be used.

Aeq = [];
beq = [];

% The lower and upper bounds of the bus power need to be set. First, it is
% desired that the historic dispatch (as recorded on SCADA) is used for the
% generator buses. To allow some flexibility for optimisation, a leeway of
% +/- 1MW is given.

LB = zeros(nbus,1);
UB = zeros(nbus,1);

for i = 1:(ngen+nhvdc)
    mask = MPCdc.bus(:,BUS_I) == MPCdc.gen(i,GEN_BUS);
    LB(mask) = LB(mask) + MPCdc.gen(i,PG) - 1;
    UB(mask) = UB(mask) + MPCdc.gen(i,PG) + 1;
end

% AEMO only releases the demand aggregated by state. To estimate the demand
% at each bus, required the statistical assumption that demand is
% proportional to population. Hence, the bus demands are approximate and to
% fine tune this estimate a 10% leeway is given. Essentially, the
% optimisation reallocates the demands so as to minimise the angle
% differences, thereby making the linearisation more valid.


for i = 1:nbus
    if MPCdc.bus(i,PD) > 0
        LB(i) = LB(i) - 1.1*MPCdc.bus(i,PD);
        UB(i) = UB(i) - 0.9*MPCdc.bus(i,PD);
    else
        LB(i) = LB(i) - 0.9*MPCdc.bus(i,PD);
        UB(i) = UB(i) - 1.1*MPCdc.bus(i,PD);        
    end
end


%% Running optimisation
Popt = quadprog(H,f,A,b,Aeq,beq,LB,UB);

%% Updating the MPC structure
% The Popt contains the optimal power injections at each bus. Now, the
% voltage angle at each bus can be calculated.

MPCdc.bus(:,VA) = pinv(Incidence'*D*Incidence)*(Popt/MPCdc.baseMVA);
MPCdc.bus(:,VA) = MPCdc.bus(:,VA) - MPCdc.bus(MPCdc.bus(:,BUS_TYPE)==REF, VA);


% Next the bus power demand and dispatch need to be updated
for i = 1:nbus
    mask = MPC.gen(:,GEN_BUS) == MPC.bus(i,BUS_I);
    MPCdc.bus(i,PD) = median([MPCdc.bus(i,PD)*[0.9,1.1], sum(MPC.gen(mask,PG)) - Popt(i)]);
    if any(mask)
        MPCdc.gen(mask,PG) = MPC.gen(mask,PG) + ...
                (Popt(i)-sum(MPC.gen(mask,PG))+MPCdc.bus(i,PD))/sum(mask);
    end
end


%%
% Pnoopt = -MPCdc.bus(:,PD);
% for i = 1:(ngen+nhvdc)
%     mask = MPC.bus(:,BUS_I) == MPC.gen(i,GEN_BUS);
%     Pnoopt(mask) = Pnoopt(mask) + MPC.gen(i,PG);
% end


%% Plotting Angles
% As a comparison, the voltage angles are plotted for the original case,
% the optisimised case and the Pyomo optimised case.

aus = shaperead("aus.shp");
load edge_geotbl
pyomoang = readtable('vang_pyomo.csv');
pyomoang = pyomoang(~nodesTAS, :);
nooptang = pinv(Incidence'*D*Incidence)*(0.5*(LB+UB)/MPCdc.baseMVA);
nooptang = nooptang - nooptang(MPCdc.bus(:,BUS_TYPE)==REF);


figure;
for i = 1:3
    subplot(1,3,i);
    plot(aus(1).X, aus(1).Y, 'k'); hold on;
    plot(aus(2).X, aus(2).Y, 'k'); hold on;
    plot(aus(3).X, aus(3).Y, 'k');
    plot(aus(4).X, aus(4).Y, 'k');
    axis('square');
    for j = 1:nline
        plot(edge_geotbl.Longitude{i}, edge_geotbl.Latitude{i}, 'Color', [0.7 0.7 0.7]);
    end
end

subplot(1,3,1);
title('Quadprog');
colormap('jet');
scatter(node_tbl.LONGITUDE, node_tbl.LATITUDE, 5, MPCdc.bus(:,VA)*180/pi, 'filled');
caxis(minmax([MPCdc.bus(:,VA); pyomoang{:,2}; nooptang]')*180/pi);

subplot(1,3,2);
title('Pyomo');
scatter(node_tbl.LONGITUDE, node_tbl.LATITUDE, 5, pyomoang{:,2}*180/pi, 'filled');
caxis(minmax([MPCdc.bus(:,VA); pyomoang{:,2}; nooptang]')*180/pi);

subplot(1,3,3);
title('No Optimisation');
scatter(node_tbl.LONGITUDE, node_tbl.LATITUDE, 5, nooptang*180/pi, 'filled');
caxis(minmax([MPCdc.bus(:,VA); pyomoang{:,2}; nooptang]')*180/pi);

colorbar('Position', [0.1, 0.25, 0.8, 0.02], 'Orientation', 'horizontal');
set(subplot(1, 3, 1), 'Position', [0.1, 0.25, 0.25, 0.7]);
set(subplot(1, 3, 2), 'Position', [0.4, 0.25, 0.25, 0.7]);
set(subplot(1, 3, 3), 'Position', [0.7, 0.25, 0.25, 0.7]);

%% Statistics on the Angles
% For better comparison, some statistics on the angles are plotted.

figure;
subplot(2,1,1);
boxplot([MPCdc.bus(:,VA), pyomoang{:,2}, nooptang]*180/pi, ...
    'Labels', {'Quadprog', 'Pyomo', 'No optimisation'});
title('Voltage Angles');

subplot(2,1,2);
boxplot(Incidence*[MPCdc.bus(:,VA), pyomoang{:,2}, nooptang]*180/pi, ...
    'Labels', {'Quadprog', 'Pyomo', 'No optimisation'});
title('Voltage Angle Differences');

%% Comparing Dispatch
% The dispatch should match the historic dispatch.

figure;
plot(MPCdc.gen(:,PG),MPC.gen(:,PG), '.'); hold on;
plot(xlim, xlim, 'r');
xlabel('Optimised dispatch');
ylabel('Historic dispatch');

%% Plotting Demand
figure;
for i = 1:3
    subplot(1,3,i);
    plot(aus(1).X, aus(1).Y, 'k'); hold on;
    plot(aus(2).X, aus(2).Y, 'k'); hold on;
    plot(aus(3).X, aus(3).Y, 'k');
    plot(aus(4).X, aus(4).Y, 'k');
    axis('square');
    for j = 1:nline
        plot(edge_geotbl.Longitude{i}, edge_geotbl.Latitude{i}, 'Color', [0.7 0.7 0.7]);
    end
end

subplot(1,3,1);
title('Quadprog');
colormap('jet');
scatter(node_tbl.LONGITUDE, node_tbl.LATITUDE, 5, MPCdc.bus(:,PD), 'filled');
colorbar;
%caxis(minmax([MPCdc.bus(:,VA); pyomoang{:,2}; nooptang]')*180/pi);

subplot(1,3,2);
title('Historic');
scatter(node_tbl.LONGITUDE, node_tbl.LATITUDE, 5, MPC.bus(:,PD), 'filled');
colorbar;

subplot(1,3,3);
title('Difference');
scatter(node_tbl.LONGITUDE, node_tbl.LATITUDE, 5, MPCdc.bus(:,PD)-MPC.bus(:,PD), 'filled');
colorbar;

%caxis(minmax([MPCdc.bus(:,VA); pyomoang{:,2}; nooptang]')*180/pi);

%colorbar('Position', [0.1, 0.25, 0.8, 0.02], 'Orientation', 'horizontal');
%set(subplot(1, 2, 1), 'Position', [0.1, 0.25, 0.25, 0.7]);
%set(subplot(1, 2, 2), 'Position', [0.4, 0.25, 0.25, 0.7]);
