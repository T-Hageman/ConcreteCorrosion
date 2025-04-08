close all
clear all
clc

addpath(genpath('../PostProcessingReqs'))

meshFile = "../../Results/mesh_0.hdf5";
dataFile = "../../Results/results_100.hdf5";
timeDataFile = "../../Results/TimeData.hdf5";
DataNamesNodes = {"H", "OH", "Fe", "FeOH", "Na", "Cl", "O2", "ePot", ...
				  "ePot", "HER", "OER", "Corrosion", ...
				  "ePot", "HER", "OER", "Corrosion"};
for i=1:8
	GroupName{i} = "/Interior/";
end
for i=9:12
	GroupName{i} = "/Bar/";
end
for i=13:16
	GroupName{i} = "/Pit/";
end

DataNamesIPs = {};

PotBased = false;

individualfigures = false;

% Load mesh
%h5disp(meshFile)

% Load data
%h5disp(dataFile)
if (individualfigures == false)
	f = figure;
	t = tiledlayout('flow');
	t.Padding = "tight";
	t.TileSpacing = "tight";
	f.Units = "centimeters";
	f.Position = [1, 0, 50, 25];
end
read = false;
while read == false;
	try
		t = h5read(dataFile,"/time")
		
		for i=1:length(DataNamesNodes)
			X{i} = h5read(meshFile,GroupName{i}+"X");
			Y{i} = h5read(meshFile,GroupName{i}+"Y");
			Z{i} = h5read(meshFile,GroupName{i}+"Z");
		
			Data = h5read(dataFile,GroupName{i}+DataNamesNodes{i});
			D{i} = Data;
		end
		for i=1:length(DataNamesIPs)
			Xip = h5read(meshFile,GroupName{i}+"Xip");
			Yip = h5read(meshFile,GroupName{i}+"Yip");
			Zip = h5read(meshFile,GroupName{i}+"Zip");
		
			DIP{i} = h5read(dataFile,GroupName{i}+DataNamesIPs{i});
		end

		%tDataTypes = h5read(timeDataFile,"/TimeDataTypes");
		%tData = h5read(timeDataFile,"/TimeData");

		read = true;
	catch ME
		pause(2)
	end
end

if (PotBased)
	for i=1:7
		D{i} = exp(D{i});
	end
end

for i=1:8
	if (individualfigures)
		figure
	else
		nexttile
	end
	PlotNodeData(X{i}, Y{i}, Z{i}, D{i});
	title(DataNamesNodes{i});
	view(0,0)
end

DataNamesNodes{9} = "E_m-ePot";
EM = tData(end,5);
D{9} = D{9} - EM;
D{13} = D{13} - EM;
for i=9:12
	if (individualfigures)
		figure
	else
		nexttile
	end
	PlotNodeData(X{i}, Y{i}, Z{i}, D{i});
	hold on
	PlotNodeData(X{i+4}, Y{i+4}, Z{i+4}, D{i+4});
	title(DataNamesNodes{i});
	view(0,90)
	axis equal
end

% if (individualfigures)
% 	figure
% else
% 	nexttile
% end

% D{1}(D{1}<0) = nan;
% pH = -log10(D{1}/1000);
% 	PlotNodeData(X, Y, Z, pH);
% 	title("pH");
% 	view(0,0)
% 
% 
% 
% if (individualfigures)
% 	figure
% else
% 	nexttile
% end
% 
% D{2}(D{2}<0) = nan;
% pOH = -log10(D{2}/1000);
% 	PlotNodeData(X, Y, Z, pOH);
% 	title("pOH");
% 	view(0,0)

% if (individualfigures)
% 	figure
% else
% 	nexttile
% end
% 
% PlotNodeData(X, Y, Z, 14-(pH+pOH));
% title("pH Error");
% view(0,0)

% if (individualfigures)
% 	figure
% else
% 	nexttile
% end
% D{1}(X>0.01) = nan;
% pH = -log10(D{1}/1000);
% 	PlotNodeData(X, Y, Z, pH);
% 	title("pH");
% 	view(0,0)

if (individualfigures)
	figure
else
	nexttile
end
yyaxis left
semilogy(tData(:,1),abs(tData(:,2)),'DisplayName','H')
hold on
semilogy(tData(:,1),abs(tData(:,3)),'DisplayName','O')
semilogy(tData(:,1),abs(tData(:,4)),'DisplayName','Fe')
xlabel('time [s]');
ylabel('I [A]');

yyaxis right
plot(tData(:,1),tData(:,5),'DisplayName','E_m')
ylabel('E_m [V_SHE]')
leg = legend('Location','southoutside','NumColumns',4);


if (individualfigures)
	figure
else
	nexttile
end
yyaxis left
plot(tData(:,1)/3600,tData(:,2),'DisplayName','H')
hold on
plot(tData(:,1)/3600,tData(:,3),'DisplayName','O')
plot(tData(:,1)/3600,tData(:,4),'DisplayName','Fe')
plot(tData(:,1)/3600,tData(:,2)+tData(:,3)+tData(:,4),'k','DisplayName','I_tot')
xlabel('time [hour]');
ylabel('I [A]');

yyaxis right
plot(tData(:,1)/3600,tData(:,5),'DisplayName','E_m')
ylabel('E_m [V_SHE]')
leg = legend('Location','southoutside','NumColumns',4);







% for i=1:length(DataNamesIPs)
% 	nexttile
% 	PlotIPData(Xip, Yip, DIP{i});
% 	title(DataNamesIPs{i});
% 	axis image
% 	if (i>1 && i<=4)
% 		clim([-2e6 0.5e6]);
% 	end
% 	if (i==5)
% 		clim([-1e6 1e6]);
% 	end
% end


function PlotNodeData(X, Y, Z, Data)
if (size(X,1)==10)
	i=0;
	for el=1:size(X,2)
		for f=1:4
			switch f
				case 1
					order = [1 2 3];
				case 2
					order = [1 2 4];
				case 3
					order = [2 3 4];
				case 4
					order = [3 1 4];
			end
			i=i+1;
			X_el(i,:) = X(order,el);
			Y_el(i,:) = Y(order,el);
			Z_el(i,:) = Z(order,el);
			Data_el(i,:) = Data(order,el);
		end
	end
	p = patch(X_el',Y_el',Z_el',Data_el','EdgeColor','interp','FaceColor','k','FaceAlpha',1.0);
	colorbar
	clim([min(min(Data_el)) max(max(Data_el))+1e-12])
elseif (size(X,1)==6)
	order = [1 2 3];
	for i=1:size(X,2)
		X_el(i,:) = X(order,i);
		Y_el(i,:) = Y(order,i);
		Z_el(i,:) = Z(order,i);
		Data_el(i,:) = Data(order,i);
	end
	p = patch(X_el',Y_el',Z_el',Data_el','EdgeColor','interp','FaceColor','interp','FaceAlpha',1.0);
	colorbar
	%clim([min(min(Data_el)) max(max(Data_el))+1e-12])
else
	if (size(X,1)==9)
		order = [1 2 3 6 9 8 7 4];
	elseif (size(X,1)==6)
		order = [1 2 3];
	else
		order = [1 2 4 3];
	end
	for i=1:size(X,2)
		X_el(i,:) = X(order,i);
		Y_el(i,:) = Y(order,i);
		Data_el(i,:) = Data(order,i);
	end
	p = patch(X_el',Y_el',Data_el',Data_el','EdgeColor','interp','FaceColor','interp');
	colorbar
	end
end

function PlotIPData(X, Y, Data)
	[xplot, yplot] = meshgrid(linspace(min(min(X)), max(max(X)),200), linspace(min(min(Y)), max(max(Y)),200));

	F = scatteredInterpolant(X(:), Y(:), Data(:));
	zplot = F(xplot, yplot);

	surf(xplot, yplot, zplot,'EdgeColor','none','FaceColor','interp');
	view(2)
	colorbar
end