%-------------------------------------
% BEGIN: function infHorizonDae.m
%-------------------------------------
function [dae Ddae]= ramseygpopsDae(sol)


t = sol.time;
k = sol.state;
c = sol.control;
p = sol.parameter;
global alpha delta n gx

%Restrição dinâmica
kdot = k.^alpha-c-(delta+n+gx).*k;

dae = (kdot);
% if nargout == 2
% dfdx = alpha.*k.^(alpha-1)-(delta+n+gx);
% dfdu = -ones(size(c));
% dfdt = zeros(size(t));
% dfdp = zeros(size(p));
% Ddae = [ dfdx dfdu dfdt dfdp];
end
% Ddae = [df/dx, df/du, df/dt, df/dp]
%-----------------------------------
% END: function infHorizonDae.m
%-----------------------------------

