function [ derivatives ] = RocketODE( 
%
%
%
%
%
%
% --------- (States In order)------------------
% 1- Mass of rocket;
% 2- Masso of Air
% 3- Volume;
% 4- Velocity x;
% 5- Velocity z;
% 6- Range (X location);
% 7- Height (Z location);
% 

% PREDEFINE Constants:


%% Phase 1: 
if 

Pressure = ( ( V0 ./ States(1) ) .^ GammaGas ) .* (Pgage+Pamb) ; 
Thrust = 2.*cd.* At .* ( PressurePhase1 - Pamb) ;
TotalVeloc = sqrt( (States(3).^2) + (States(4).^2) );
Drag = ( RohAirBoulder / 2) .* (TotalVeloc).^2 * CD*Ab; 


%% Phase 2:
elseif i=2
    
    
%% Phase 3: 

else
    
end


end