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

global t2 t1 t3
%globals are where the time of each phase ends, so last elment in t1 is
%where phase 1 ends, so on.

%% assign golbals:

t1 = [0 0];
t2 = [0];
t3 = [0];
%% 

%initial conditions:

VelX0 = 0;
VelZ0 = 0;
Range0 = 0;
Height0 = y0;

% Call ODE
[ Time Results ] = ode45(@(Time,States) RocketODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,y0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R), [ 0:0.001:5],[TotalMass0 MassAirInit...
VAirInit VelX0 VelZ0 Range0 y0 ]);



%% calculate thrust

%ODE Has really a hard way of calculating the time, even though matrices t1
%and t2 can tell us the end time of each phase there's no time that matches
%it exactly in the Time matrix, and extracting the force of thrust via
%globals really yields bad results, thus the forces has to be calculated
%manually.


%index namings: those will use to index the storing of thrust and time
store1 = 0;
store2 = 0;
store3 = 0;

[ r c ] = size(Results);
for i=1:r
    
if Results(i,3) < Volbottle
    
Pressure1(i) = ( ( VAirInit ./ Results(i,3) ) .^ GammaGas ) .* (Pgage+Pamb) ; 
Thrust1(i) = 2.* Cd .* ThroatArea .* ( Pressure1(i) - Pamb) ;
TP1(i) = Time(i);

%phase 2
elseif Results(i,3)>= Volbottle
    
    %T and P of end states
Tend = TAirInit * (( VAirInit/Volbottle) ^ (GammaGas-1) );
Pend = (Pgage+Pamb) * (( VAirInit/Volbottle) ^ (GammaGas) );
    
PressureCond = Pend * (Results(i,2)/MassAirInit)^(GammaGas) ;

if PressureCond>Pamb
Density = Results(i,2) / Volbottle;
Temp = PressureCond/(Density*R);
CriticalP = (PressureCond) * (2./(GammaGas+1)).^(GammaGas/(GammaGas-1));
  
if CriticalP > Pamb
    
    Mach  = 1;
    Texit = (2/(GammaGas+1))*Temp ;
    Vexit = sqrt(GammaGas*Texit*R);
    Pexit = CriticalP;
    Densityexit = CriticalP/(R*Texit) ;
    
elseif CriticalP <= Pamb
    
   Mach = sqrt(( (PressureCond/Pamb)^( ( (GammaGas-1)/GammaGas)) - 1 ) * (2/(GammaGas-1)));
   Texit = Temp/(1+((GammaGas-1)/2)*Mach^2);
   Pexit = Pamb;
   Densityexit = Pamb/(R*Texit) ;
   Vexit = Mach * sqrt(GammaGas*Texit*R);
end

% how mass of air and rocket change with time

MassAirFlowRate = Cd*Densityexit*ThroatArea*Vexit;

store2 = store2 + 1;
Thrust2(store2) = MassAirFlowRate *Vexit + (Pexit-Pamb)*ThroatArea ;
TP2(store2) = Time(i);

else
store3 = store3 + 1;
Thrust3(store3) = 0 ;
TP3(store3) = Time(i);


%% Phase 3: 
end




end

end


%% plot thrust:

figure(1);
plot([TP1 TP2 TP3],[Thrust1 Thrust2 Thrust3],'Color',[1 0.5 0.2],'LineWidth',1.4)
hold on
plot(TP1(end),Thrust1(end),'*','Color',[ 0 0.5 0],'MarkerSize',7,'MarkerFaceColor',[0 0.5 0])
plot(TP2(end),Thrust2(end),'*','Color',[0.2 0 0],'MarkerSize',7,'MarkerFaceColor',[0.2 0 0])
xlim([0 TP3(90)])
grid minor
title('Thrust VS Time')
xlabel('Time (Seconds)')
ylabel('Thrust (N)')
legend('Thrust profile','End of Water Phase','End of Air Phase','Location','NorthEast')


%% for future devlopments, ignore.

%{
% this portion of the graph was taken from MIT matlab training section,

% http://web.mit.edu/8.13/matlab/MatlabTraining_IAP_2012/AGV/DemoFiles/ScriptFiles/html/Part3_Animation.html

% figure;
% 
% subplot(2,1,1);
% xlabel('time (sec)'); ylabel('angle (\circ)');
% 
% for id = 1:length(Time3)
%    subplot(2,1,1);
%    plot(Results(id,6),Results(id,7), 'LineWidth', 2);
%    line(Results(id,6), Results(id,7), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
%    line(Results2(id,6), Results2(id,7), 'Marker', '*', 'MarkerSize', 20, 'Color', 'r');
%    line(Results3(id,6), Results3(id,7), 'Marker', '^', 'MarkerSize', 20, 'Color', 'b');
% 
% %    line(Time(id), Results(id,7), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
%    xlabel('Range (m)'); ylabel('Height (m)');
% 
% %    % The bottom plot shows the animation of the double pendulum
% %    subplot(2,1,2);
% %    plot([0, x(id,1);x(id,1), x(id,2)], [0, y(id,1);y(id,1), y(id,2)], ...
% %       '.-', 'MarkerSize', 20, 'LineWidth', 2);
%    axis equal; axis([0 max(Results(:,6))+2 0 max(Results(:,7))+2]);
%    title(sprintf('Time: %0.2f sec', Time(id)));
%    grid minor
% 
%    drawnow;
%    
% end

%}
%{
clf;

ax = axes('XLim',[0 max(Results(:,6))+2],'YLim',[ 0 max(Results(:,7))+2],...
    'ZLim',[ 0 1 ]);

view(3)
grid minor;
axis equal

% create the body of the rocket:

% [ xcone ycone zcone ] = cylinder([0.2 -2]); % cylinder with top raduis of 0 so it's cone
% [ xcylin ycylin zcylin ] = cylinder([0.2 0.05]);

[ xcone ycone zcone ] = cylinder([0.4 0.05]); % cylinder with top raduis of 0 so it's cone
[ xcylin ycylin zcylin ] = cylinder([0.4 0.4]);

%create body parts of the rocket, assumption: no fins, we don't need them $_$

Parts(1) = surface(xcone,ycone,zcone,'FaceColor','red');
Parts(2) = surface(xcylin,ycylin,1.2*zcylin-1.2,'FaceColor','yellow');


%put the rocket together!
%SET BREAKPOINT HERE! ERROR IN THE NEXT LINE!!
Rocket = hgtransform('Parent',ax);
set(Parts,'Parent',Rocket)

%}