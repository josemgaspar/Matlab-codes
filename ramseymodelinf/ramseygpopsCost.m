function [ Mayer, Lagrange, DMayer, DLagrange ] = ramseygpopsCost( sol )
%Funcional objetivo

global rho sigma gx n



t0 = sol.initial.time;
k0 = sol.initial.state;
tf = sol.terminal.time;
kf = sol.terminal.state;
t = sol.time;
k = sol.state;
c = sol.control;
p = sol.parameter;

Mayer = zeros(size(t0));
%Integrando
Lagrange = - (exp(-(rho-n).*t)).*(((c.*exp(gx.*t)).^(1-(1/sigma))-1)/(1-(1/sigma)));

% if nargout == 4
% % DMayer = [ dM/dx0, dM/dt0, dM/dxf,
% DMayer = [zeros(1,length(k0)), zeros(1,length(t0)), zeros(1,length(kf)), ...
% ... % dM/dtf, dM/dp]
% zeros(1,length(tf)), zeros(1,length(p))];
% % DLagrange = [ dL/dx, dL/du, dL/dp, dL/dt]
% DLagrange =[ zeros(size(k), - (exp(-(rho-n).*t)).*exp(gx.*t).^(1-(1/sigma)).*...
%     (1-(1/sigma)).*(c.^(-(1/sigma)))/(1-(1/sigma))), zeros(length(t),length(p)),...
%     (exp(t.*(n - rho)).*(n - rho).*((c.*exp(gx.*t)).^(1 - 1/sigma) - 1))...
%     /(1/sigma - 1) - (c.*gx.*exp(gx.*t).*exp(t.*(n - rho)))./(c.*exp(gx.*t))^(1/sigma)];
% end
end

