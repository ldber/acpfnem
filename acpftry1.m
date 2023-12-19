%% Flat start with historic data
runpf(MPC); %Does not converge!!!


%% Initialise angles from DC power flow results
MPC.bus(:,VA) = nooptang;
runpf(MPC); %Does not converge!!!

runpf(MPCdc); %Does not converge!!!

%% Set reactive power demands for 1pu voltages at each bus
Ybus = (Incidence.')*diag(1./(MPC.branch(:,BR_R) + 1j*MPC.branch(:,BR_X)))*Incidence + ...
    diag(sum(diag(1j*MPC.branch(:,BR_B)/2)*abs(Incidence)));

MPC.bus(:,VA) = 0;
Vbus = MPC.bus(:,VM).*exp(1j*MPC.bus(:,VA));
Sbus = Vbus.*(conj(Ybus*Vbus));
MPC.bus(:,QD) = -imag(Sbus);
runpf(MPC); %Does not converge!!!

MPC.bus(:,VA) = nooptang;
Vbus = MPC.bus(:,VM).*exp(1j*MPC.bus(:,VA));
Sbus = Vbus.*(conj(Ybus*Vbus));
MPC.bus(:,QD) = -imag(Sbus);
runpf(MPC); %Does not converge!!!

Vbus = MPCdc.bus(:,VM).*exp(1j*MPCdc.bus(:,VA));
Sbus = Vbus.*(conj(Ybus*Vbus));
MPCdc.bus(:,QD) = -imag(Sbus);
runpf(MPCdc); %Does not converge!!!

% Ybus = (Incidence.')*diag(1./((edge_tbl.R_PU + 1j*edge_tbl.X_PU)./edge_tbl.NUM_LINES))*Incidence + ...
%     diag(sum(diag(1j*(edge_tbl.B_PU).*edge_tbl.NUM_LINES/2)*abs(Incidence)));
% Pbus = zeros(nbus,1);
% Qbus = zeros(nbus,1);
% for i = 1:nbus
%     for k = 1:nbus
%         Pbus(i) = Pbus(i) + MPC.bus(i,VM)*MPC.bus(k,VM)*abs(Ybus(i,k))*cos(MPC.bus(i,VA)-MPC.bus(k,VA)-angle(Ybus(i,k)));
%         Qbus(i) = Qbus(i) + MPC.bus(i,VM)*MPC.bus(k,VM)*abs(Ybus(i,k))*sin(MPC.bus(i,VA)-MPC.bus(k,VA)-angle(Ybus(i,k)));
%     end
% end

%% Set line resistance to zero
MPCr0 = MPCdc;
MPCr0.branch(:,BR_R) = 0;
MPCr0.bus(:,QD) = 0;
MPCr0.bus(:,VA) = 0;
MPCr0 = runpf(MPCr0); %Converged!!!

%% Plotting Voltage Magnitudes and Angle Differences
aus = shaperead("aus.shp");
figure;
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
title('Voltage Magnitudes');
colormap('jet');
scatter(node_tbl.LONGITUDE, node_tbl.LATITUDE, 5, MPCr0.bus(:,VM), 'filled');
colorbar;

subplot(1,2,2);
title('Voltage Angle Difference');
colormap('jet');
scatter(node_tbl.LONGITUDE, node_tbl.LATITUDE, 5, (MPCr0.bus(:,VA)-MPCdc.bus(:,VA)), 'filled');
colorbar;

%% Plotting Power Flows

