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

TempSample = (T_Sample_1+T_Sample_2)/2;


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

maxIndi = find(max(TempSample)==TempSample);
[ m3 b3 sig_y3 sig_b3 sig_m3 Q3 ] = LSM(Time(maxIndi:end),TempSample(maxIndi:end));

TimeSampleAdded = Time(235);


%[ m3 b3 sig_y3 sig_b3 sig_m3 Q3 ] = LSM(Time(300:end),TempSample(300:end));
%[ m3 b3 sig_y3 sig_b3 sig_m3 Q3 ] = LSM(Time(280:end),TempSample(280:end));


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

Temp1 = mean(T_boiling(1:235)); % the temp T1
Temp1_unc = std(T_boiling(1:235))/sqrt(length(T_boiling(1:235))); %Uncertinity with it.

%% Uncertainty measurements

% Get uncertinnty with new values along the fitting line

% establish the matrices that will have each new uncertinnty as we step
% away from our best fit, this's just for us so we can see what's going on.

sigma_newY1 = zeros(1,length(Time));
sigma_newY2 = zeros(1,length(Time));
sigma_newY3 = zeros(1,length(Time));


for i=1:length(Time)
    
    sigma_newY1(i) = [ Time(i) 1 ] * Q1 * [ Time(i) ; 1 ];
    sigma_newY2(i) = [ Time(i) 1 ] * Q2 * [ Time(i) ; 1 ];
    sigma_newY3(i) = [ Time(i) 1 ] * Q3 * [ Time(i) ; 1 ];
    
end


%uncertininty in T2 (Equilibrium temp)

sigmaT2 = [ TimeT2 1 ] * Q1 * [ TimeT2 ; 1 ];

sigmaT2 = double(sqrt(sigmaT2));

%% Specific Heat

SpecificHeatSample = (SpecifHeatCalo*Calo_mass*(Temp2-Temp_L)) / ((Sample_mass*(Temp1-Temp2)));

%convert units

SpecificHeatSample = SpecificHeatSample * ( 1 /0.238846 );

%% error estimations 

%-=-=-=-==-=-=-=-=-=-=-=-=-=-=(General Method)=-=-=-=-=-=-=-=-=-=


%-=-=-=-==-=-=-=-=-=-=-=-=-=-=(Step by step)=-=-=-=-=-=-=-=-=-=

% A = T2 - T0 (T0 == Temp_L)

A = Temp2 - Temp_L ;
sigmaA = ( ( sigmaT2 ) ^2 + ( sig_y1 ) ^2 ) ^(1/2);


% B = T1 - T2
B = (Temp1 - Temp2);
sigmaB =  ( ( sigmaT2 ) ^2 + ( Temp1_unc ) ^2 ) ^(1/2);

% A/B = D

D = A/B;
sigmaD = abs(D) * ( ( sigmaB/B ) ^2 + ( sigmaA/A ) ^2 ) ^(1/2);

% C = mc / ms ( mass of calorimeter / mass of sample )

C = Calo_mass/Sample_mass ;
sigmaC = abs(C) * ( ( unc_Calo_mass/Calo_mass ) ^2 + ( unc_Sample_mass/Sample_mass ) ^2 ) ^(1/2);

% E = C * D

E = C * D;
sigmaE = abs(E) *  ( ( sigmaC/C ) ^2 + ( sigmaD/D ) ^2 ) ^(1/2);

% final uncertinty, E * Specific Heat of calorimeter since it's treated as
% exact.

SigmaSpecificHeat = sigmaE * (SpecifHeatCalo*( 1 /0.238846 ));
SigmaSpecificHeat = double(SigmaSpecificHeat);

%% print out the results:

fprintf('Initial temperature of calorimeter is: %f \n',Temp_L);
fprintf('Initial water temp when sample was added: %f \n',T_boiling(235));
fprintf('Time when the sample was added (seconds) is: %f \n',TimeSampleAdded);
fprintf('Halfway Temp is: %f \n',Temp_mid);
fprintf('Equilibrium temp of the sample and calorimete is: %f \n',Temp2);
fprintf('\n');
fprintf('--------------------------------------------------------- \n');
fprintf('Your Sample Specific Heat is: %f \n',SpecificHeatSample);
fprintf('with uncertainty of: %f',SigmaSpecificHeat);
fprintf(',which is: %f',(SigmaSpecificHeat/SpecificHeatSample)*100);
fprintf('%% %f \n');
fprintf('');

%% plot the analysis results

figure(1)

scatter(Time,TempSample,2,'*','MarkerEdgeColor',[0.7 0.9 0.6])
hold on
plot(Time,output_line_fit1,'--b','LineWidth',1)
hold on
plot(Time,output_line_fit2,'--r','LineWidth',1)
hold on
plot(Time,output_line_fit3,'--k','LineWidth',1)
hold on
plot([TimeSampleAdded TimeSampleAdded], [0 40],'b','LineWidth',1.2)
hold on
plot([TimeT2 TimeT2], [0 40],'k','LineWidth',1.2)
hold on
plot(TimeSampleAdded,Temp_L,'r<','LineWidth',3)
hold on
plot(TimeT2,Temp_mid,'ro','LineWidth',3)
hold on
plot(TimeT2,Temp2,'r*','LineWidth',3)
hold on
grid minor
ylim([20 28])
legend('Collected Data','Best fit 1','Best fit 2','Best fit 3','V-line when sample was added','V-line halfway','Sample Initial temperature','Midway sample temperature','Equilibrium temperature')
title('Temperature Profile For Calorimeter')
xlabel('Time (Seconds)')
ylabel('Temperature (celsius)')

%% plot the uncertininty with possible samples.


figure(2)

%diffrence boxes

box_x_1=[mean(Pb) SpecificHeatSample SpecificHeatSample mean(Pb)];
box_y_1=[0.4 0.4 1.3 1.3];

box_x_2=[ SpecificHeatSample Tellurium_Copper Tellurium_Copper SpecificHeatSample ];
box_y_2=[0.4 0.4 1.3 1.3];


box_x_3=[ Tellurium_Copper Zn_Cu_Ti Zn_Cu_Ti Tellurium_Copper ];
box_y_3=[0.4 0.4 1.3 1.3];

box_x_4=[  Zn_Cu_Ti Al_6063_T1 Al_6063_T1 Zn_Cu_Ti  ];
box_y_4=[0.4 0.4 1.3 1.3];

%plot
patch(box_x_1,box_y_1,'red','FaceAlpha',0.08)
hold on
patch(box_x_2,box_y_2,'green','FaceAlpha',0.08)
hold on
patch(box_x_3,box_y_3,'yellow','FaceAlpha',0.08)
hold on
patch(box_x_4,box_y_4,'blue','FaceAlpha',0.08)
hold on
herrorbar(Zn_Cu_Ti,1,0.00001,0.00001,'^')
hold on
herrorbar(mean(Pb),1,(std(Pb)/length(Pb)), (std(Pb)/length(Pb)),'o')
hold on
herrorbar(Tellurium_Copper,1,0.00001,0.00001,'*')
hold on
herrorbar(Al_6063_T1,1,0.00001,0.00001,'h')
hold on
herrorbar(SpecificHeatSample,0.8,SigmaSpecificHeat,SigmaSpecificHeat,'p')
title('Candidate samples')
ylim([0.4 1.3])
xlabel('Specific Heat with uncertainties')
set(gca, 'YTick', []) %hide y-axis
grid minor
legend('Difference between Pb & Sample', 'Difference between Sample & Copper','Difference between Copper & Zn-Cu-Ti ','Difference between Zn-Cu-Ti & Al', '','Zn-Cu-Ti ','','Pb','','Cu','','Al','','Sample')

