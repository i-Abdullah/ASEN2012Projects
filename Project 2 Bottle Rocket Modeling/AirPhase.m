function Derivatives2 = AirPhase(time,States2,Tend,Pend,Voltbottle,R,GammaAir,Pamb,ThroatArea,RohAirBoulder,MassAirInit,cd,CD,Ab,g)
%
% ------------------------ ( STATES IN ORDER ) ------------------
%   1- Mass of Air
%   2- Mass (of the whole rocket), it changes as the water is expelled
%   3- Velocity x
%   4- Velocity z
%   5- Position x
%   6- Position z
%   7- Pressure
%
%
%

Pressure2 = Pend * (MassAirInit/States2(1))^(GammaAir) ;
Density2 = States2(1) / Voltbottle;
Temp2 = Pressure2/(Density2*R);
CriticalP = (Pressure2) * (2./(GammaAir+1)).^(GammaAir/(GammaAir-1));

if CriticalP > Pamb
    Mach  = 1;
    Texit = (2/(GammaAir+1))*Temp2 ;
    Vexit = sqrt(GammaAir*Texit*R);
    Pexit = CriticalP;
    Densityexit = CriticalP/(R*Texit) ;
    
elseif CriticalP <= Pamb
    
   Mach = sqrt(((Pressure2/Pamb)^(( (GammaAir-1)/GammaAir) - 1 )) * (2/(GammaAir-1)))
   Texit = Temp2*(1+((GammaAir-1)/2)*Mach^2)
   Pexit = Pamb;
   Densityexit = Pamb/(R*Texit) ;
   Vexit = Mach * sqrt(GammaAir*Texit*R);
    
end

% define everything again :(

TotalVeloc = sqrt( (States2(3).^2) + (States2(4).^2) );
HeadingX = ( (States2(3)) / TotalVeloc );
HeadingZ =  ( (States2(4)) / TotalVeloc );
States2(7) = Pexit;


%% Derivatives

%derivative of mass: how massdecreases
MassAirFlowRate = cd*Densityexit*ThroatArea*Vexit;
MassRocketFlowRate = -MassAirFlowRate;

% thrust and drag and total force

Drag2 = ( RohAirBoulder / 2) .* (TotalVeloc).^2 * CD*Ab ; 
Thrust2 = MassAirFlowRate *Vexit + (Pexit-Pamb)*ThroatArea ;

%Accerlation (Velocity Derivative)
DAccelration_Dt_InX = ( (( Thrust2 - Drag2) * HeadingX)) ./ States2(2) ;
DAccelration_Dt_InZ =  ( ((Thrust2 - Drag2) * HeadingZ) - States2(2)*g ) ./ States2(2) ;


%%

Derivatives2 = [ MassAirFlowRate ; MassRocketFlowRate ; DAccelration_Dt_InX ; DAccelration_Dt_InZ; States2(3) ; States2(4) ; 0 ]
   