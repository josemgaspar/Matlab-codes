function [ ydot ] = tobinblanchard( t,y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global sigma a1 a2 g k u b1 m0
ydot=zeros(2,1);
ydot(1)= sigma*(a1-1)*y(1) + sigma*a2*y(2)+sigma*g;
ydot(2)= ((k*y(2))/u-b1)*y(1)-(y(2)*m0)/u;
end