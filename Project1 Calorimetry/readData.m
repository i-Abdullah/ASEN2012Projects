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

[ m1 b1 ] = LSM(Time(1:235),T_Sample(1:235));
[ m2 b2 ] = LSM(Time(235:280),T_Sample(235:280));
%option 1 [ m3 b3 ] = LSM(Time(340:end),T_Sample(340:end));
%option 2 [ m3 b3 ] = LSM(Time(280:end),T_Sample(280:end));
%option 3
[ m3 b3 ] = LSM(Time(300:end),T_Sample(300:end));
TimeSampleAdded = Time(235);

%-=-=-=-=-=-=-=-=-=-=-=(Find T1)=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



coeff1 = [ m1 ; b1 ] ;
coeff2 = [ m2 ; b2 ] ;
coeff3 = [ m3 ; b3 ] ;
output_line_fit1 = polyval(coeff1,Time);
output_line_fit2 = polyval(coeff2,Time);
output_line_fit3 = polyval(coeff3,Time);

f1 = @(x) m1*x +b1;
f2 = @(x) m2*x +b2;
f3 = @(x) m3*x +b3;

%% get temp values
%vertical line @ the time the sample was added in
x = Time(235);

T_Equation1 = @(x) f1(x) - f2(x);


TimeTemp1 = fzero(T_Equation1,0);
TimeTemp2 = fzero(T_Equation2,0);
TimeTemp3 = (TimeTemp1 + TimeTemp2)/2;

Temp_1 = feval(f1,TimeTemp1);
Temp_2 = feval(f3,TimeTemp2);
Temp_3 = (Temp_1 + Temp_2)/2;


fprintf('Initial temperature of calorimeter is: %f \n',Temp_1);
fprintf('Time when the sample was added (seconds) is: %f \n',TimeTemp2);
fprintf('Equilibrium temp of the sample and calorimete is: %f \n',Temp_3);
%% plot
scatter(Time,T_Sample)
hold on
plot(Time,output_line_fit1,'--','LineWidth',2)
hold on
plot(Time,output_line_fit2,'--','LineWidth',2)
hold on
plot(Time,output_line_fit3,'--r','LineWidth',2)
hold on
line([TimeTemp2 TimeTemp2], [0 40])
hold on
plot(Time(235),T_Sample(235),'r*')
hold on
plot(TimeTemp1,Temp_1,'r*')
hold on
plot(TimeTemp2,Temp_2,'r*')
hold on
plot(TimeTemp3,Temp_3,'r*')
hold on
grid minor
ylim([20 30])