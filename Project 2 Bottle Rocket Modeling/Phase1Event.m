function [value, isterminal, direction] = Phase1Event(time,States)
%this function is used to set-up stoping criteria for ode45.
value = States(1) - 0.002 ;
isterminal = 1;   % Stop the integration
direction  = 0;

end