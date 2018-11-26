%% info


%% houskeeping

clear;
clc;
close all;


%% inital conditions 

global NetForcePhase1

g = 9.81; % m/s2, acceleration due to gravity,
Cd= 0.8; % discharge coefficient
Rohairamb = 0.961; % kg/m^3 ambient air density
Volbottle= 0.002; % m^3 volume of empty bottle
Pamb= 12.1*6894.76; % converted to Pa, atmospheric pressure
GammaAir = 1.4; % ratio of specific heats for air
RohWater = 1000; % kg/m^3, density of water
DThroat= 2.1; % cm, diameter of throat
DBottle= 10.5; % in cm, diameter of bottle
R = 287; %J/kgK, gas constant of air
MBottle= 0.15; % kg mass of empty 2-liter bottle with cone and fins
CD= 0.5; % drag coefficient
Pgage= 50*6894.76; % in pascal, the 6894.76 is to convert. initial gage pressure of air in bottleVolwater,
VWaterInit= 0.001; % m^3, initial volume of water inside bottle
TAirInit = 300; % K, initial temperature of
Airv0 = 0.0 ;% m/s, initial velocity of rocket
Theta= pi/4 ; % initial angle of rocket
X0 = 0.0; % in meters, initial horizontal distance
y0 = 0.25; % in m, initial vertical height
TestStandLength= 0.5; % in m, length of test stand
VAirInit = Volbottle - VWaterInit ; %initial volume of Air.
ThroatArea = pi * ((DThroat*10^-2)/2)^2;
BottleArea  = pi * ((DBottle*10^-2)/2)^2;
TotalMassInit = MBottle + (VWaterInit*RohWater) + (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit));
RohAirBoulder = 0.961;
NetForcePhase1 = []; %First column is thrust, second is drag, third is pressure inside bottle
%% Numerical integration.

%% PHASE 1: WATER

%initial conditions:

VelXInit = 0;
VelZInit = 0;
RangeInit = 0;
HeightInit = 0;


% event stoping condition will be used, it's in seperate function.
Stop = odeset('Events',@Phase1Event);

[ Time States1,TimeEventHappens1,ValueOfEvent1 ] = ode45(@(time,States) WaterPhase(time,States,Cd,...
ThroatArea,RohWater,(Pgage+Pamb),VAirInit,Pamb,GammaAir,Pgage,CD,BottleArea,...
RohAirBoulder,g), [ 0 20],[VAirInit TotalMassInit...
 VelXInit VelZInit RangeInit HeightInit ],Stop);

%% PHSE 2: Air

%Initial Values:
[ r c ] = size(States1) ;
[ r1 c1 ] = size(Time);

MassRocket2Init = States1(r,2);
MassAir2Init = MassAirInit;
Velcoity2_X_Init = States1(r,3);
Velcoity2_Z_Init = States1(r,4);
Position2_X_Init = States1(r,5);
Position2_Z_Init = States1(r,6);
Pressure2Init = NetForcePhase1(r1,3);
TimeStart =  Time(r1,1);

Stop2 = odeset('Events',@Phase2Event);

Tend = TAirInit * (( VAirInit/Volbottle) ^ (GammaAir-1) );
Pend = (Pgage+Pamb) * (( VAirInit/Volbottle) ^ (GammaAir) );

[ Time States2,TimeEventHappens2,ValueOfEvent2 ] = ode45 (@(time,States2) AirPhase(time,States2,Tend,Pend,Volbottle,R,GammaAir,Pamb,ThroatArea,RohAirBoulder,MassAirInit,Cd,CD,BottleArea,g)...
    ,[TimeStart 20],[ MassAir2Init;MassRocket2Init;Velcoity2_X_Init;Velcoity2_Z_Init;Position2_X_Init;Position2_Z_Init;Pressure2Init],Stop2);