MPCq0 = MPCr0;

J = makeJac(ext2int(MPCq0), 1);
dQ_dVrms = full(J(nbus+1:2*nbus,nbus+1:2*nbus));

Qbus = -100*dQ_dVrms*(1-MPCq0.bus(:,VM));

Qbus1 = median([-abs(MPCq0.bus(:,PD))/0.9*sqrt(1-0.9^2),abs(MPCq0.bus(:,PD))/0.9*sqrt(1-0.9^2), Qbus], 2);

% for i = 1:nbus
%     if MPCq0.bus(i,VM)>1
%         Qbus(i) = min(MPCq0.bus(i,VM)-1, 0.05) * (abs(MPCq0.bus(i,PD))/0.9*sqrt(1-0.9^2))/0.05;
%     else
%         Qbus(i) = -min(1-MPCq0.bus(i,VM), 0.05) * (abs(MPCq0.bus(i,PD))/0.9*sqrt(1-0.9^2))/0.05;        
%     end
% end

%% Plotting voltages
figure;
for i = 1:2
    subplot(1,2,i)
    plot(aus(1).X, aus(1).Y, 'k'); hold on;
    plot(aus(2).X, aus(2).Y, 'k');
    plot(aus(3).X, aus(3).Y, 'k');
    plot(aus(4).X, aus(4).Y, 'k');
    axis('square');
    for k = 1:nline
        plot(edge_geotbl.Longitude{k}, edge_geotbl.Latitude{k}, 'Color', [0.7 0.7 0.7]);
    end
end

subplot(1,2,1)
title('Voltage RMS');
colormap('jet');
scatter(node_tbl.LONGITUDE, node_tbl.LATITUDE, 5, MPCq0.bus(:,VM), 'filled');
colorbar;

subplot(1,2,2)
title('Reactive Power Demand');
colormap('jet');
scatter(node_tbl.LONGITUDE, node_tbl.LATITUDE, 5, Qbus, 'filled');
colorbar;

%%
MPCq0.bus(:,QD) = Qbus1;
MPCq0 = runpf(MPCq0);

figure;
plot(aus(1).X, aus(1).Y, 'k'); hold on;
plot(aus(2).X, aus(2).Y, 'k');
plot(aus(3).X, aus(3).Y, 'k');
plot(aus(4).X, aus(4).Y, 'k');
axis('square');
for k = 1:nline
    plot(edge_geotbl.Longitude{k}, edge_geotbl.Latitude{k}, 'Color', [0.7 0.7 0.7]);
end

title('New Voltage RMS');
colormap('jet');
scatter(node_tbl.LONGITUDE, node_tbl.LATITUDE, 5, MPCq0.bus(:,VM), 'filled');
colorbar;

