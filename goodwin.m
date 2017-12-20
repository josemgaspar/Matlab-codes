function [ ydot ] = goodwin( t,y )
%Defini��o de vari�veis globais
global alpha beta delta v RHO
%ydot � o vector das solu��es do modelo goodwin:
%y(1) � a "share" ou peso dos sal�rios dos funcion�rios na economia;
%y(2) � a taxa de emprego na economia;
%delta � a taxa de crescimento da produtividade;
%RHO � o retorno m�dio de capital para todos os investimentos
%v d�-nos a varia��o da produtividade ao longo do tempo sobre a popula��o
%alhpa e beta s�o par�metros tais que a taxa de crescimento dos sal�rios �
%igual a -alhpa*beta*y(1)*y(2)
ydot=[-(delta+alpha)*y(1)+beta*y(1)*y(2); 
      (1/RHO-delta-v)*y(2)-(y(1)*y(2))/RHO];
end

