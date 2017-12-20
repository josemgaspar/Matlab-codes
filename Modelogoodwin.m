clc;
close all; clear all;
global alpha beta delta v RHO
alpha=0.4;beta=0.9;delta=0.05;v=0.01;RHO=5;
y1i=0.5; y2i=0.5;time=[0,100];
[t,ydot]=ode23s(@goodwin,time,[y1i y2i]);
plot(t,ydot); title('Comportamento económico cíclico - Modelo Goodwin');
xlabel('tempo');ylabel('"Share" dos funcionários na economia; Taxa de emprego');
legend('share dos funcionários','taxa de emprego');