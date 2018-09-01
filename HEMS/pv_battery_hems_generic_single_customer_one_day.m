%% *********************************************************************************************************************************************************
% ELEC5206 Sustainable Energy Systems Group Assignment: Techno-economic feasibility study of a small-scale PV-battery system
% Mixed integer linear program for a residential PV-battery system
% Donald Azuatalam (donald.azuatalam@sydney.edu.au)
% Dr. Gregor Verbic and Dr. Archie Chapman
%% *********************************************************************************************************************************************************

% clear
% clear class
clear cplex
close all
clc


ampl = AMPL; % Create an AMPL instance
ampl.reset;
ampl.setOption('solver', 'cplex');   % set solver
ampl.eval('option version;')    % Display AMPL version
basef = fileparts(which('pv_battery_hems_generic_single_customer_one_day'));  % Find matlab file
ampl.read([basef  '\' 'pv_battery_hems.mod']);  % Find ampl file

N = 48; %Number of time slots in a day
Cust = 1; %Number of customers
Days = 1; % One day simulation
etaI = 1; %Inverter efficiency
etaBc = sqrt(0.84); etaBd = etaBc; %Battery charge and discharge efficiency
D = (1:96)';  % Half-hours in 2 days (2-day rolling horizon)
D = num2cell(D);
days = ampl.getSet('D');
days.setValues(D);

%PV_Size = [5 4 5 3 6 6 5 4 4 9 3 4 8 4 6 5 10 5 4 4 5 4 7 4 4 4 5 7 7 5 5 4 5 3 5 5 5 10 5 5 7 5 5 3 6 3 4 5 5 5 4 4]; %PV Sizes for 52 customers
G1 = [6.5 9.8 14]; % Battery maximum SOC for LG Chem 6.5, LG Chem 10, Tesla Power Wall II
G2 = [0.6 1.0 0.5]; % Battery minimum SOC
G3 = [4.2 5.0 5.0]; % Battery maximum charge and discharge rates
G4 = sqrt([0.84 0.84 0.84]); % Battery charge and discharge efficiencies
G5 = G1/2; % Start of day 1 SOC

%GENERIC DATA
input_file = 'Data.xlsx';
%Demand data
PD = xlsread(input_file,'','B2:B49'); % Electrical demand (daily load profile)
PD((N*Days)+1:(N*Days)+N) = PD((N*Days)-(N-1):(N*Days)); %dummy profile for second day (used only for rolling horizon)
PD(PD<0) = 0;
%PV data
PV = xlsread(input_file,'','C2:C49'); % PV output (daily PV output profile)
PV((N*Days)+1:(N*Days)+N) = PV((N*Days)-(N-1):(N*Days));
PV(PV<0) = 0;
PV_Size = round(max(PV)); %approximate PV size of generic PV profile
% Time of Use Tariff (Ausgrid) - Retail electricity tariff
c_g = xlsread(input_file,'','D2:D49');
c_g =  repmat(c_g,2,1);
c_g = num2cell(c_g);
% Modify feed-in-Tariff from AMPL model file

%Pre-allocate memory for results
Pgplus = zeros(N*Days,Cust);
Pgminus = zeros(N*Days,Cust);
Pbplus = zeros(N*Days,Cust);
Pbminus = zeros(N*Days,Cust);
sb = zeros(N*Days,Cust);
dg = zeros(N*Days,Cust);
eb = zeros(N*Days,Cust);


%Assigning battery specs to PV Size

%Check PV size of customer and assign a corresponding battery size
if (PV_Size <= 4) % Use LG Chem 6.5
    ampl.getParameter('ebM').setValues(G1(1));
    ampl.getParameter('ebm').setValues(G2(1));
    ampl.getParameter('PbM').setValues(G3(1));
    ampl.getParameter('etaBc').setValues(G4(1))
    ampl.getParameter('eb1').setValues(G5(1));
    ebM = G1(1);ebm = G2(1);PbM = G2(1);
elseif (PV_Size > 4) && (PV_Size <= 6) % Use LG Chem 10
    ampl.getParameter('ebM').setValues(G1(2));
    ampl.getParameter('ebm').setValues(G2(2));
    ampl.getParameter('PbM').setValues(G3(2));
    ampl.getParameter('etaBc').setValues(G4(2))
    ampl.getParameter('eb1').setValues(G5(2));
    ebM = G1(2);ebm = G2(2);PbM = G2(2);
elseif (PV_Size > 6) % Use Tesla Power Wall II
    ampl.getParameter('ebM').setValues(G1(3));
    ampl.getParameter('ebm').setValues(G2(3));
    ampl.getParameter('PbM').setValues(G3(3));
    ampl.getParameter('etaBc').setValues(G4(3))
    ampl.getParameter('eb1').setValues(G5(3));
    ebM = G1(3);ebm = G2(3);PbM = G2(3);
end

for i = 1:Days  %Solve for one day
    
    disp(i);
    
    U = PD;
    V = PV;
    
    U = U(N*(i-1)+1:N*(i+1)); % Solve 2 days at a time (2 day rolling horizon)
    V = V(N*(i-1)+1:N*(i+1));
    
    df = DataFrame(1, 'D', 'Pd');  % Create a dataframe with 1 index for Demand
    df.setMatrix(U, D);
    ampl.setData(df);
    
    df = DataFrame(1, 'D', 'Ppv');  % Create a dataframe with 1 index for PV
    df.setMatrix(V, D);
    ampl.setData(df);
    
    %Necessary only for time-of-use tariffs
    df = DataFrame(1, 'D', 'c_g'); % % Create a dataframe with 1 index for TOU Tariff
    df.setMatrix(c_g, D);
    ampl.setData(df);
    
    ampl.eval('option presolve 1;'); %Enable presolve, otherwise 0
    ampl.eval('option cplex_options ''mipgap=1e-1'';');
    
    ampl.solve() %Solve
    
    cost = ampl.getObjective('cost'); % Get objective map by AMPL name
    %fprintf('Objective is: %f\n', cost.value());  % Print it
    
    % GET VARIABLES
    varPgplus = ampl.getVariable('Pgplus').getValues('val').getColumnAsDoubles('val'); %Get values of Pgplus in a dataframe object %Power import from grid
    Pgplus(N*(i-1)+1:N*i) = varPgplus(1:N); % For 2 day rolling horizon optimization, only keep results for day 1
    
    varPgminus = ampl.getVariable('Pgminus').getValues('val').getColumnAsDoubles('val');  % Get values of Pgminus in a dataframe object %Power export to grid
    Pgminus(N*(i-1)+1:N*i) = varPgminus(1:N);
    
    varPbplus = ampl.getVariable('Pbplus').getValues('val').getColumnAsDoubles('val'); %Get values of Pgplus in a dataframe object %Battery charge power
    Pbplus(N*(i-1)+1:N*i) = varPbplus(1:N);
    
    varPbminus = ampl.getVariable('Pbminus').getValues('val').getColumnAsDoubles('val');  % Get values of Pgminus in a dataframe object %Battery discharge power
    Pbminus(N*(i-1)+1:N*i) = varPbminus(1:N);
    
    varsb = ampl.getVariable('sb').getValues('val').getColumnAsDoubles('val');  % Get the values of eb in a dataframe object %Battery state of charge
    sb(N*(i-1)+1:N*i) = varsb(1:N);
    
    vardg = ampl.getVariable('dg').getValues('val').getColumnAsDoubles('val');  % Get the values of eb in a dataframe object %Battery state of charge
    dg(N*(i-1)+1:N*i) = vardg(1:N);
    
    vareb = ampl.getVariable('eb').getValues('val').getColumnAsDoubles('val');  % Get the values of eb in a dataframe object %Battery state of charge
    eb(N*(i-1)+1:N*i) = vareb(1:N);
    ampl.getParameter('eb1').setValues(vareb(N)); %End of day 1 (of 2) SOC == Beginning of next 2 day's SOC (Rolling horizon approach)
    
end


%Save results in csv files
csvwrite('MILP_generic_one_day_variable_Pgplus.csv',Pgplus);
csvwrite('MILP_generic_one_day_variable_Pgminus.csv',Pgminus);
csvwrite('MILP_generic_one_day_variable_Pbplus.csv',Pbplus);
csvwrite('MILP_generic_one_day_variable_Pbminus.csv',Pbminus);
csvwrite('MILP_generic_one_day_variable_sb.csv',sb);
csvwrite('MILP_generic_one_day_variable_dg.csv',dg);
csvwrite('MILP_generic_one_day_variable_eb.csv',eb);
Pgnet = Pgplus-Pgminus; %Net Grid Power


%Visualise results
n = 1:1*N;
NN = 1*N;
PV = PV(n);
PD = PD(n);
%% N/B: You may need to adjust axes limits (ax.XLim and ax.YLim) if the need arises

fig1 = figure;
fig1.Units = 'centimeters';
fig1.Position = [10   10   16   11.1125];
fig1.Color = 'white'; % set(gcf,'color','white')

subplot(3,2,1)
stairs(n,Pgnet)
title('Grid power')
ax = gca; ax.FontSize = 8; ax.XTickLabel = {'0','6','12','18','24'}; ax.XTick = 2*[0 6 12 18 24];
ax.XLim = [1 NN]; ax.YLim = [0 5]; ax.Box = 'on'; ax.Color = 'white';
xlabel('Hours')
ylabel('{\itP} (kW)')

subplot(3,2,2)
plot(n,eb);
title('Baterry state of charge')
ax = gca; ax.FontSize = 8; ax.XTickLabel = {'0','6','12','18','24'}; ax.XTick = 2*[0 6 12 18 24];
ax.XLim = [1 NN]; ax.YLim = [0 7]; ax.Box = 'on'; ax.Color = 'white';
xlabel('Hours')
ylabel('{\itE} (kWh)')

subplot(3,2,3); hold on
stairs(n,Pbplus);
stairs(n,Pbminus);
title('Battery charge/discharge power')
ax = gca; ax.FontSize = 8; ax.XTickLabel = {'0','6','12','18','24'}; ax.XTick = 2*[0 6 12 18 24];
ax.XLim = [1 NN]; ax.YLim = [0 5]; ax.Box = 'on'; ax.Color = 'white';
xlabel('Hours')
ylabel('{\itP} (kW)')
legend('charge','discharge')

subplot(3,2,4)
plot(n,PD)
title('Electrical demand')
ax = gca; ax.FontSize = 8; ax.XTickLabel = {'0','6','12','18','24'}; ax.XTick = 2*[0 6 12 18 24];
ax.XLim = [1 NN]; ax.YLim = [0 2]; ax.Box = 'on'; ax.Color = 'white';
xlabel('Hours')
ylabel('{\itP} (kW)')

subplot(3,2,5)
plot(n,PV)
title('PV output')
ax = gca; ax.FontSize = 8; ax.XTickLabel = {'0','6','12','18','24'}; ax.XTick = 2*[0 6 12 18 24];
ax.XLim = [1 NN]; ax.YLim = [0 1.6]; ax.Box = 'on'; ax.Color = 'white';
xlabel('Hours')
ylabel('{\itP} (kW)')

c_g = xlsread(input_file,'','D2:D49')'; % Retail electricity tariff
subplot(3,2,6)
stairs(c_g)
title('Retail tariff')
ax = gca; ax.FontSize = 8; ax.XTickLabel = {'0','6','12','18','24'}; ax.XTick = 2*[0 6 12 18 24];
ax.XLim = [1 NN]; ax.YLim = [0 0.5]; ax.Box = 'on'; ax.Color = 'white';
xlabel('Hours')
ylabel('{\itc} ($/kWh)')

saveas(fig1, 'house_variables.fig');
export_fig house_variables.pdf

%% 
fig2 = figure;
fig2.Units = 'centimeters';
fig2.Position = [5   5   16   16];
fig2.Color = 'white'; % set(gcf,'color','white')

subplot(4,2,1)
stairs(n,Pgplus)
title('Pgplus')
ax = gca; ax.FontSize = 8; ax.XTickLabel = {'0','6','12','18','24'}; ax.XTick = 2*[0 6 12 18 24];
ax.XLim = [1 NN]; ax.YLim = [-0.1 5]; ax.Box = 'on'; ax.Color = 'white';
xlabel('Hours')
ylabel('{\itP} (kW)')

subplot(4,2,2)
stairs(n,Pgminus)
title('Pgminus')
ax = gca; ax.FontSize = 8; ax.XTickLabel = {'0','6','12','18','24'}; ax.XTick = 2*[0 6 12 18 24];
ax.XLim = [1 NN]; ax.YLim = [-0.1 5]; ax.Box = 'on'; ax.Color = 'white';
xlabel('Hours')
ylabel('{\itP} (kW)')

subplot(4,2,3)
stairs(n,Pbplus)
title('Pbplus')
ax = gca; ax.FontSize = 8; ax.XTickLabel = {'0','6','12','18','24'}; ax.XTick = 2*[0 6 12 18 24];
ax.XLim = [1 NN]; ax.YLim = [-0.1 5]; ax.Box = 'on'; ax.Color = 'white';
xlabel('Hours')
ylabel('{\itP} (kW)')

subplot(4,2,4)
stairs(n,Pbminus)
title('Pbminus')
ax = gca; ax.FontSize = 8; ax.XTickLabel = {'0','6','12','18','24'}; ax.XTick = 2*[0 6 12 18 24];
ax.XLim = [1 NN]; ax.YLim = [-0.1 5]; ax.Box = 'on'; ax.Color = 'white';
xlabel('Hours')
ylabel('{\itP} (kW)')

subplot(4,2,5)
stairs(n,eb)
title('eb')
ax = gca; ax.FontSize = 8; ax.XTickLabel = {'0','6','12','18','24'}; ax.XTick = 2*[0 6 12 18 24];
ax.XLim = [1 NN]; ax.YLim = [-0.1 14]; ax.Box = 'on'; ax.Color = 'white';
xlabel('Hours')
ylabel('{\itE} (kWh)')

subplot(4,2,6)
stairs(n,sb)
title('sb')
ax = gca; ax.FontSize = 8; ax.XTickLabel = {'0','6','12','18','24'}; ax.XTick = 2*[0 6 12 18 24];
ax.XLim = [1 NN]; ax.YLim = [-0.1 1.1]; ax.Box = 'on'; ax.Color = 'white';
xlabel('Hours')
ylabel('0/1')

subplot(4,2,7)
stairs(n,dg)
title('dg')
ax = gca; ax.FontSize = 8; ax.XTickLabel = {'0','6','12','18','24'}; ax.XTick = 2*[0 6 12 18 24];
ax.XLim = [1 NN]; ax.YLim = [-0.1 1.1]; ax.Box = 'on'; ax.Color = 'white';
xlabel('Hours')
ylabel('0/1')


saveas(fig2, 'decision_variables.fig');
export_fig decision_variables.pdf
%%

fig3 = figure; hold on;
fig3.Units = 'centimeters';
fig3.Position = [10   10   16   8];
fig3.Color = 'white'; % set(gcf,'color','white')

ax = gca; % current axes
ax.FontSize = 8;
ax.XTickLabel = {'0','6','12','18','24'};
ax.XTick = 2*[0 6 12 18 24];
ax.XLim = [1 NN];
ax.YLim = [0 2];
ax.Box = 'on';
ax.Color = 'white'; % set(gca,'color','white')


%%Visualising Results
%%
[from_pv_demand,from_pv_battery,from_pv_grid,from_grid_battery,from_battery_demand,from_grid_demand,from_battery_grid] = extra_variables(NN,PV,PD,Pgplus,Pgminus,Pbplus,Pbminus,etaI,etaBc,etaBd);
PV_check = (from_pv_demand + from_pv_battery + from_pv_grid)/etaI; 
PD_check = from_grid_demand + from_battery_demand + from_pv_demand;
total_demand = [from_grid_demand,from_battery_demand,from_pv_demand,from_pv_battery,from_pv_grid];
area_stairs(n', total_demand)


stairs(n,etaI*PV)
stairs(n,from_pv_demand)
title('Household electrical demand and PV generation')
xlabel('Hours')
ylabel('{\itP} (kW)')
source = {'From Grid to Demand','From Battery to Demand','From PV to Demand','From PV to Battery','From PV to Grid','PV','PV to Demand'};
legend(source, 'Location', 'NorthWest')
colormap summer
%%

grid on
grid minor

saveas(fig3, 'house_demand.fig');
export_fig house_demand.pdf