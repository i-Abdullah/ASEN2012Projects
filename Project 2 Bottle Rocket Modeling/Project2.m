%% info


%% houskeeping

clear;
clc;
close all;


%% inital conditions 


g = 9.81; % m/s2, acceleration due to gravity,
Cd= 0.8; % discharge coefficient
Rhoairamb = 0.961; % kg/m^3 ambient air density
Volbottle= 0.002; % m^3 volume of empty bottle
Pamb= 12.1*6894.76; % converted to Pa, atmospheric pressure
GammaGas = 1.4; % ratio of specific heats for air
RhoWater = 1000; % kg/m^3, density of water
DThroat= 2.1; % cm, diameter of throat
DBottle= 10.5; % in cm, diameter of bottle
R = 287; %J/kgK, gas constant of air
MBottle= 0.15; % kg mass of empty 2-liter bottle with cone and fins
CD= 0.5; % drag coefficient
Pgage= 50*6894.76; % in pascal, the 6894.76 is to convert. initial gage pressure of air in bottleVolwater,
VWaterInit= 0.001; % m^3, initial volume of water inside bottle
TAirInit = 300; % K, initial temperature of
Airv0 = 0.0 ;% m/s, initial velocity of rocket
Theta= 45 ; % initial angle of rocket in degress
X0 = 0.0; % in meters, initial horizontal distance
y0 = 0.25; % in m, initial vertical height
TestStandLength= 0.5; % in m, length of test stand
VAirInit = Volbottle - VWaterInit ; %initial volume of Air.
ThroatArea = pi * ((DThroat*10^-2)/2)^2; %Area of throat
BottleArea  = pi * ((DBottle*10^-2)/2)^2; %Bottle Area
TotalMass0 = MBottle + (VWaterInit*RhoWater) + (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit));

%% Numerical integration.

%% 

%initial conditions:

VelX0 = 0;
VelZ0 = 0;
Range0 = 0;
Height0 = y0;

[ Time Results ] = ode45(@(Time,States) RocketODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,y0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R), [ 0 5],[TotalMass0 MassAirInit...
VAirInit VelX0 VelZ0 Range0 y0 ]);

%% plot thrust

Pressure = ( ( VAirInit ./ Results(:,3) ) .^ GammaGas ) .* (Pgage+Pamb) ;
Thrust = 2.* Cd .* ThroatArea .* ( Pressure - Pamb) ;

subplot(2,1,1)
comet(Results(:,6),Results(:,7))
ylim([0 40])
subplot(2,1,2)
comet(Results(:,6),Results(:,7))
ylim([0 40])

grid minor
