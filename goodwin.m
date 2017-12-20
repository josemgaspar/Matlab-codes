function [ ydot ] = goodwin( t,y )
%Definição de variáveis globais
global alpha beta delta v RHO
%ydot é o vector das soluções do modelo goodwin:
%y(1) é a "share" ou peso dos salários dos funcionários na economia;
%y(2) é a taxa de emprego na economia;
%delta é a taxa de crescimento da produtividade;
%RHO é o retorno médio de capital para todos os investimentos
%v dá-nos a variação da produtividade ao longo do tempo sobre a população
%alhpa e beta são parâmetros tais que a taxa de crescimento dos salários é
%igual a -alhpa*beta*y(1)*y(2)
ydot=[-(delta+alpha)*y(1)+beta*y(1)*y(2); 
      (1/RHO-delta-v)*y(2)-(y(1)*y(2))/RHO];
end

