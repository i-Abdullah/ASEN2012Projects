%% read the file

clear
clc
close all;

Data = load('Sample_A.txt'); % load the file
Time = Data(:,1); % time
T_boiling = Data(:,2); % boiling temp


T_Sample_1 = Data(:,3); % temp of sample using thermocouple 1
T_Sample_2 = Data(:,4); % temp of sample using thermocouple 2

%-=-=-=-=-=-=-=-=-=-=-=-=-= ( Avg temp between 1 and 2 )%-=-=-=-=-=-=-=-=-=-=-=-=

TempSample = (T_Sample_1+T_Sample_2)/2


%-=-=-=-=-=-=-=-=-=-=-=-=-= ( Material info )%-=-=-=-=-=-=-=-=-=-=-=-=

% from the file / Manually

Sample_mass = 91.767; %in grams
unc_Sample_mass = 0.001; %uncertainty

Calo_mass = 318.3; %in grams
unc_Calo_mass = 0.05;  %uncertainty

SpecifHeatCalo = 0.214;
%-=-=-=-=-=-=-=-=-=-=-=-=-= ( Possible Materials )%-=-=-=-=-=-=-=-=-=-=-=-=
%specific heats that are given so we can compare for different alloys.

Zn_Cu_Ti = 0.402;
Tellurium_Copper = 0.261;
Pb = 0.100386:0.001:0.129;
Al_6063_T1 = 0.9;



%% Linear Fitting to find T0, Temp when the 

[ m1 b1 sig_y1 sig_b1 sig_m1 Q1 ] = LSM(Time(1:235),TempSample(1:235));
[ m2 b2 sig_y2 sig_b2 sig_m2 Q2 ] = LSM(Time(235:280),TempSample(235:280));
[ m3 b3 sig_y3 sig_b3 sig_m3 Q3 ] = LSM(Time(300:end),TempSample(300:end));
TimeSampleAdded = Time(235);


%option 1 [ m3 b3 sig_y3 sig_b3 sig_m3 ] = LSM(Time(340:end),TempSample(340:end));
%option 2 [ m3 b3 sig_y3 sig_b3 sig_m3 ] = LSM(Time(280:end),TempSample(280:end));


%-=-=-=-=-=-=-=-=-=-=-=(Establish the fit lines)=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

% get the coefficients
coeff1 = [ m1 ; b1 ] ;
coeff2 = [ m2 ; b2 ] ;
coeff3 = [ m3 ; b3 ] ;

%evaluate them along the time interval to establish a line.
output_line_fit1 = polyval(coeff1,Time); %from t=0 to t=235 seconds
output_line_fit2 = polyval(coeff2,Time); %from t=235 to t=280 seconds
output_line_fit3 = polyval(coeff3,Time); %from t=280 to t=end seconds

%put them into a matlab function
f1 = @(x) m1*x +b1;
f2 = @(x) m2*x +b2;
f3 = @(x) m3*x +b3;

%% get temp values TH, TL, and Mid Temp
Temp_L = feval(f1,TimeSampleAdded);
Temp_H = feval(f3,TimeSampleAdded);
Temp_mid = (Temp_L+Temp_H)/2;

% Now get T2 where a line that passes through Temp_Mid intercepts Output_line_fit3
% this line is @ what time the temp reached temp_mid, wo we will solve for
% it using root finding methods

T2_poly = @(x) m2*x + b2 - Temp_mid; %convert t2 into root finding problem
TimeT2 = solve(T2_poly); %find the root.
Temp2 = f3(TimeT2); %the root is the time, thus evaluate that time at the third best fit
Temp2 = double(Temp2); %convert it from symbolic expression to double number.

% T1 will be the temp of of the boiling water @ the time the sample was
% inputted, thus we will take the mean from time = 1 to the time we think
% the sample was inputted and use it, and consider it as our T1, and the
% standard deviation will be the uncertintity with it.

T1 = mean(T_boiling(1:235); % the temp T1
T1_unc = mean(T_boiling(1:235); %Uncertinity with it.

fprintf('Initial temperature of calorimeter is: %f \n',Temp_L);
fprintf('Time when the sample was added (seconds) is: %f \n',TimeSampleAdded);
fprintf('Equilibrium temp of the sample and calorimete is: %f \n',Temp2);
fprintf('Halfway Temp is: %f \n',Temp_mid);
fprintf('Initial water temp when sample was added: %f \n',T_boiling(235));

%% Uncertainty measurements

% Get uncertinnty with new values along the fitting line

% establish the matrices

sigma_newY1 = zeros(1,length(Time));
sigma_newY2 = zeros(1,length(Time));
sigma_newY3 = zeros(1,length(Time));

for i=1:length(Time)
    
    sigma_newY1(i) = [ Time(i) 1 ] * Q1 * [ Time(i) ; 1 ];
    sigma_newY2(i) = [ Time(i) 1 ] * Q2 * [ Time(i) ; 1 ];
    sigma_newY3(i) = [ Time(i) 1 ] * Q3 * [ Time(i) ; 1 ];
    
end





%% Specific Heat
SpecificHeatSample = (SpecifHeatCalo*Calo_mass*(Temp2-Temp_L)) / ((Sample_mass*(mean(T_boiling(1:235))-Temp2)));
%convert units

SpecificHeatSample = SpecificHeatSample * ( 1 /0.238846 );

%% error estimations 

%% preapre errors for plot

sigmay1 = ones(1,length(TempSample))*sig_y1;
%% plot
scatter(Time,TempSample,2,'*','MarkerEdgeColor',[0.7 0.9 0.6])
hold on
%shadedErrorBar(Time(1:235),TempSample(1:235),(ones(235,1)*sig_y1))
fill([Time(1:end);flipud(Time(1:end))],[output_line_fit1(1:end)-sigma_newY1(1:end);flipud(output_line_fit1(1:end)+sigma_newY1(1:end))],[.9 .9 .9],'linestyle','none');
hold on
plot(Time,output_line_fit1,'--r','LineWidth',1)
hold on
plot(Time,output_line_fit2,'-.r','LineWidth',1)
hold on
plot(Time,output_line_fit3,'-.r')
hold on
plot([TimeSampleAdded TimeSampleAdded], [0 40],'-.b')
hold on
plot([TimeT2 TimeT2], [0 40],'-.b')
hold on
plot(Time(235),T_Sample_1(235),'r*')
hold on
plot(TimeSampleAdded,Temp_L,'b*')
hold on
plot(TimeT2,Temp_mid,'b*')
hold on
plot(TimeT2,Temp2,'b*')
hold on
grid minor
ylim([20 28])
