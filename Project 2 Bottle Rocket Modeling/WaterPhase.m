function Derivatives = WaterPhase(time,States,cd, At, RohH2O, P0, V0, Pamb, Gamma,Pgage,CD,Ab,RohAirBoulder,g,DragPhase1,NetForcePhase1)
%{ The following function is part of modeling a water bottle rocket, for the
% firs phase, where the rocket ejects water, Done by:
%
% - Brendan Palmer
% - Abdulla Al Ameri
%
% ------------------------ ( STATES IN ORDER ) ------------------
%
%   1- Volume
%   2- Mass (of the whole rocket), it changes as the water is expelled
%   3- Velocity x
%   4- Velocity z
%   5- Position x
%   6- Position z
%
%

%% Phase1

global NetForcePhase1

PressurePhase1 = ( ( V0 ./ States(1) ) .^ Gamma ) .* (Pgage+Pamb) ; 
ThrustPhase1 = 2.*cd.* At .* ( PressurePhase1 - Pamb) ;
TotalVeloc = sqrt( (States(3).^2) + (States(4).^2) );
DragPhase1 = ( RohAirBoulder / 2) .* (TotalVeloc).^2 * CD*Ab; 

%Heading vectors:


%if we're still on the stand our heading will always be defined by the sign
%45 
if sqrt(States(5).^2+(States(6)).^2) <= 0.5
    
HeadingX = cosd(45);
HeadingZ = sind(45);
    
else
    
HeadingX = ( (States(3)) / TotalVeloc );
HeadingZ =  ( (States(4)) / TotalVeloc );

end

%Derivatives:

%Volume (how volume changes with time)
DVolume_Dt = cd * At * sqrt ( (2/RohH2O) * (( P0 * ( V0/States(1) ) ^ (Gamma)) - Pamb ));

%Mass (how mass of water changes with time)
DMass_Dt = - cd .* At .* sqrt ( 2.*RohH2O.* ( PressurePhase1 - Pamb ) );

%Accerlation (Velocity Derivative)
DAccelration_Dt_InX = ( (( ThrustPhase1 - DragPhase1) * HeadingX)) ./ States(2) ;
DAccelration_Dt_InZ =  ( ((ThrustPhase1 - DragPhase1) * HeadingZ) - States(2)*g ) ./ States(2) ;

% Velocity (Velocity is position Derivatives)

% DP/DV : From equation 6
%DP_DV = (RohH2O/(4*TotalVeloc)) * ( 2*States(2) + 2*States(3) ) ;

NetForcePhase1 = [NetForcePhase1; ThrustPhase1 DragPhase1 PressurePhase1];
Derivatives = [ DVolume_Dt; DMass_Dt; DAccelration_Dt_InX ; DAccelration_Dt_InZ ; States(3) ; States(4) ] ;
end
