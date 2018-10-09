%% read the file

clear
clc
close all;

Data = load('Sample_A.txt');
Time = Data(:,1);
T_boiling = Data(:,2);
T_Sample = Data(:,3);
T_Sample_2 = Data(:,4);

% from the file/ Manually
Sample_mass = 91.767; %in grams
unc_Sample_mass = 0.001;

Calo_mass = 318.3;
ucn_Calo_mass = 0.05;

%% Linear Fitting to find T0, Temp when the 

[ m b ] = LSM_FirstDeg(Time(1:236),T_Sample(1:236));
coeff = [ m ; b ] ;
output_line_fit = polyval(coeff,Time);
scatter(Time,T_Sample)
hold on
plot(Time,output_line_fit)
hold on
plot(Time(236),T_Sample(236),'r*')
grid minor