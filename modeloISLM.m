%% Modelo IS-LM din�mico
%Fun��o que realiza uma an�lise ao modelo IS-LM (economia fechada ou grande
%economia aberta) din�mico. Incorpora a utiliza��o de um solver ode45 para
%sistemas de equa��es diferenciais para representar as solu��es do produto
%e da taxa de juro ao longo do tempo (Problema de valor inicial); representa as solu��es de y(t)
%contra as de r(t) em diagramas de fase. Para os diagramas de fase, cria
%ainda um campo vectorial com vectores de velocidade que permite inferir e
%tornar intuitivas as direc��es das traject�rias nos diagramas de fase. Representa ainda
%a evolu��o das principais vari�veis agregadas macroecon�micas e faz uma
%an�lise � estabilidade dos pontos de quil�brio do sistema e � inclina��o das curvas IS e LM.
%O output gerado s�o valores concretos para algumas das vari�veis
%mencionadas e 5 gr�ficos, sendo os primeiros 4 agrupados numa figura e o
%outro numa 2�.

function modeloISLM
clc; clear all; close all;      %Fechar todas as figuras abertas e limpar a janela de comandos
global alpha C I G X Q b T jay g q h beta M P L k u
disp('------------------------');
disp('Din�mica do modelo IS-LM'),
disp('------------------------');
%%Defini��o de vari�veis
%Flexibilizar a box gerada pelo comando inputdlg
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

%% A atribui��o de valores para as vari�veis � feita atrav�s de uma box de di�logo.

%Vari�veis que n�o afectam estabilidade do equil�brio:
prompt = {'Componente aut�noma do consumo C=','Componente aut�noma do Investimento I=',...
    'Componente aut�noma dos gastos p�blicos G=','Exporta��es X=','Componente aut�noma das importa��es Q=',...
    'Oferta nominal de moeda M=','N�vel geral de pre�os P=','Procura aut�noma de moeda L='};
name = 'Par�metros do modelo IS-LM (n�o afectam estabilidade do equil�brio)';
numlines = 1;
defaultanswer = {'6','25','15','6','9','30','1','3'};     %Valores atribu�dos por defeito - Eq� est�vel.
answer1 = inputdlg(prompt,name,numlines,defaultanswer,options);       %Faz aparecer a box para inserir os valores.

if isempty(answer1)     %Sai do programa se n�o forem atribu�dos valores (escolher cancelar).
    clc;
    return
end

C=str2double(answer1{1});       %Consumo aut�nomo
I=str2double(answer1{2});       %Investimento aut�nomo
G=str2double(answer1{3});       %Despesa P�blica aut�noma
X=str2double(answer1{4});       %Exporta��es
Q=str2double(answer1{5});       %Importa��o aut�noma
M=str2double(answer1{6});       %Oferta nominal de moeda
P=str2double(answer1{7});       %N�vel geral de pre�os
L=str2double(answer1{8});       %Procura aut�noma por moeda

%Vari�veis que afectam a estabilidade do sistema
prompt={'Par�metro da equa��o diferencial dy/dt=alpha(A+a1*y(t)-hr(t)  alpha=',...
    'Par�metro da equa��o diferencial dr/dt=beta(-M/P+L+ky(t)-ur(t))  beta=','Taxa marginal de imposto T=',...
    'Propens�o marginal a consumir b=','Sensibilidade do investimento ao produto real jay=',...
    'Sensibilidade do investimento � taxa de juro nominal h=','Sensibilidade da procura real de moeda face ao produto real k=',...
    'Sensibilidade da procura real de moeda face � taxa de juro nominal u=','Propens�o marginal � despesa p�blica g=',...
    'Propens�o marginal � importa��o q='};
name='Par�metros do modelo IS-LM (afectam a estabilidade do equil�brio)';
numlines=1;
defaultanswer={'0.5','0.4','0.25','0.7','0.3','1.3','0.25','0.3','0.05','0.1'};        %Valores atribu�dos por defeito - Eq� est�vel.
answer2=inputdlg(prompt,name,numlines,defaultanswer,options);       % ""

if isempty(answer2)     %""
    clc;
    return
end

alpha=str2double(answer2{1});       %Par�metro alpha da equa��o diferencial dy/dt
beta=str2double(answer2{2});        %Par�metro beta da equa��o diferencial dr/dt
T=str2double(answer2{3});           %Taxa de imposto marginal
b=str2double(answer2{4});           %Propens�o marginal ao consumo
jay=str2double(answer2{5});         %Sensibilidade do investimento ao produto real
h=str2double(answer2{6});           %Sensibilidade do investimento � taxa de juro nominal
k=str2double(answer2{7});           %Sensibilidade da procura real de moeda face ao produto real
u=str2double(answer2{8});           %Sensibilidade da procura real de moeda face � taxa de juro nominal
g=str2double(answer2{9});           %Propens�o marginal � despesa p�blica
q=str2double(answer2{10});          %Propens�o marginal � importa��o

%Valor inicial e final para o tempo e condi��es iniciais y(0) e r(0)
prompt={'Valor inicial do produto real  y(0)=','Valor inicial da taxa de juro nominal  r(0)=',...
    'tempo inicial  t0=','tempo final  tf='};
name='PVI: Condi��es iniciais';
numlines=1;
defaultanswer={'75','8','0','100'};        %Valores por defeito
answer3=inputdlg(prompt,name,numlines,defaultanswer,options);       % ""

if isempty(answer3)     % ""
    clc;
    return
end

y0=str2double(answer3{1});          %Valor inicial do produto (em t0)
r0=str2double(answer3{2});          %Valor inicial da taxa de juro (em t0)
t0=str2double(answer3{3});          %t inicial
tf=str2double(answer3{4});          %t final

a1=b*(1-T)+jay+g-q-1;       %Inverso do sim�trico do multiplicador keynesiano a1
A=C+G+I+X-Q;                %Despesa aut�noma
%%%

%% Apresenta��o do modelo e par�metros usados
 %Mostrar na janela de comandos os valores escolhidos para os par�metros
fprintf('\nModelo din�mico IS-LM:\n                           dy/dt=alpha(A+a1*y(t)-hr(t)\n                           dr/dt=beta(-M/P+L+ky(t)-ur(t))\n');
fprintf('\nCurvas IS e LM:\n                    (IS)r=(A+a1*y)/h\n                    (LM)r=(k*y-M/P+L)/u\n ');
fprintf('\ny: Produto real\nr: Taxa de juro nominal\nA: Despesa aut�noma = C+G+I+X-Q = %g\na1: b*(1-T)+jay+g-q-1 = %g\n',...
    A,a1);
fprintf('\nPar�metros usados:\nC = %g; I =  %g; G = %g; X = %g; Q = %g; M = %g; P = %g; L = %g; alpha = %g; beta = %g;\n',...
    C,I,G,X,Q,M,P,L,alpha,beta);
fprintf('T = %g; b = %g; jay = %g; h = %g; k = %g; u = %g; g = %g; q = %g\n------------------------------------------------------------\n',...
    T,b,jay,h,k,u,g,q);                                                    

%O programa d� erro se o produto de equil�brio for negativo
if k*h-a1*u<0
    errordlg('Produto de equil�brio � negativo','Erro da Fun��o');
end

%% PVI: Problema de valor inicial
[t,ydot]=ode45(@islmdinamico,[t0 tf],[y0 r0]);      %O solver ode45 atribuir a ydot os 
%valores das solu��es das equa��es diferenciais definidas em
%@islmdinamico(ver no final).

%Produto e/ou taxa de juro n�o devem ser negativos em nenhum momento do
%tempo
if all(ydot(:) > 0)        %Se todos os valores do produto e da taxa de juro forem positivas, continua.
else
    %Cria uma box em que � dado a escolher continuar com a fun��o ou n�o,
    %caso y ou r negativos em algum tempo
    ButtonName = questdlg('O valor escolhido para os par�metros implica n�veis do produto e/ou taxa de juro negativos em algum momento do tempo. Continuar simula��o?', ...
                         'Aten��o!','Sim', 'N�o','N�o');
                     fprintf('\n!!!Produto e/ou taxa de juro n�o podem ser inferiores ou iguais a zero em nenhum momento do tempo:\nSimula��o meramente ilustrativa da din�mica do modelo. \n\n')
                     if strcmp(ButtonName,'N�o')==1     %Cancela o programa
                         fprintf('\nSimula��o cancelada: Produto e/ou taxa de juro n�o podem ser inferiores ou iguais a zero em nenhum momento do tempo\n\n');
                             return
                    end
end
%%%


%% C�lculo do ponto de equil�brio do sistema de equa��es diferenciais(steady-state).
options=optimset('Display','off');      %Suprimir informa��es do solver fsolve
[pe]=fsolve(@(y) islmdinamico(t,y),[y0;r0],options);        %Calcula y e r de equil�brio (solu��o local) dado dy/dt= e dr/dt=0
%e utilizando y0 e r0 como condi��es iniciais.
disp(' ');
disp('Steady-State:');
fprintf('\nProduto de equil�brio    Taxa de juro de equil�brio\n      %f                 %f\n',pe(1,1),pe(2,1));        %Retorna o ponto de eq�

%
%% Estudo da estabilidade do ponto de equil�brio
% Estudo da estabilidade em fun��o dos valores pr�prios da matriz jacobiana
Jac=[alpha.*a1 -alpha*h;beta*k -beta*u];        %Matriz jacobiana no ponto de equil�brio.
disp(' ');
disp('Valores pr�prios da matriz jacobiana');
disp(eig(Jac));        %Mostra os val pp da matriz jacobiana


if real(eig(Jac))<0 
    if isreal(eig(Jac))==0  %Para ser est�vel, ambos os valores pr�prios t�m que ter parte real negativa
    disp('O equil�brio � um foco est�vel.');
    else
        disp('O equil�brio � um n� est�vel.')
    end
elseif real(eig(Jac))>0 
    if isreal(eig(Jac))==0
    disp('O equil�brio � um foco inst�vel.')
    elseif isreal(eig(Jac))==1
    disp('O equil�brio � um n� inst�vel.')
    end
elseif real(eig(Jac))==0 
    if isreal(eig(Jac))==0
    disp('O equil�brio � um ponto centro est�vel. Existem solu��es peri�dicas');
    end
else    
    disp('O equil�brio � inst�vel.')        %N�o � possivel haver pontos sela, pois det(Jac)>0, j� que queremos um y de equil�brio positivo.
end

if a1<0
    fprintf('\nEstabilidade do equil�brio devido a:\nMultiplicador keynesiano maior que zero (a1<0)\n');
elseif a1>0 && alpha*a1<beta*u      %Rela��o entre a1 e Multiplicador e consequ�ncias dos seus sinais (e/ou do tra�o da jacobiana) sobre a estabilidade.
    fprintf('\nEstabilidade do equil�brio devido a:\nMultiplicador keynesiano menor que zero (a1>0) e tr(Jac) < 0\n(nota:a1=b*(1-T)+jay+g-q-1\n');
else 
    fprintf('\nInstabilidade do equil�brio devido a:\nMultiplicador keynesiano menor que zero (a1>0) e tr(Jac) > 0\n(nota:a1=b*(1-T)+jay+g-q-1\n');
end

fprintf('\nValor de outras vari�veis agregadas no steady-state:\n');     %Organiza��o e reprodu��o d output na janela de comando

%Consumo, Investimento, Despesa p�blica, Exporta��es e Importa��es no
%steady-state: Output na janela de comandos
Ceq=C+b*(1-T)*pe(1,1);
Ieq=I-h*pe(2,1)+jay*pe(1,1);
Imp=Q+q*pe(1,1);
Geq=G+g*pe(1,1);
fprintf('\nConsumo    Investimento    Despesa P�blica     Exporta��es      Importa��es\n %.3f        %.3f          %.3f             %.3f            %.3f\n',...
    Ceq,Ieq,Geq,X,Imp);
disp('--------------------------------------------------------------------------------');
%%%

%% Multiplicadores do modelo IS-LM est�tico
mult1=-1/a1;     %Multiplicador keynesiano
mult2=-1/(a1-h*k/u);        %Multiplicador dos gastos p�blicos
mult3=-(k/u)/(a1-h*k/u);      %Multiplicador monet�rio
mult4=b*pe(1,1)/(a1-h*k/u);      %Multiplicador da taxa de imposto
fprintf('\nMultiplicadores (Modelo IS-LM est�tico):\n \nMultiplicador keynesiano (-1/a1):  %f\nMultiplicador dos gastos p�blicos: %f \nMultiplicador monet�rio:           %f \nMultiplicador da taxa de imposto: %f\n',...
    mult1,mult2,mult3,mult4);
disp('--------------------------------------------');
disp(' ');
disp('Inclina��o das curvas IS e LM:');     
disp(' ');
if a1<0         %A inclina��o pode ser definida em fun��o do multiplicador keynesiano (rela��o com a1)
    disp('(Multiplicador keynesiano maior que zero (a1<0): IS com inclina��o negativa)');
elseif a1==0
    disp('(Multiplicador keynesiano indeterminado (a1=0): IS com inclina��o nula; IS horizontal)');
else
    disp('(Multiplicador keynesiano menor que zero (a1>0): IS com inclina��o positiva)');
end
fprintf('\nInclina��o IS (a1/h)    Inclina��o LM (k/u)\n      %f              %f\n',a1/h,k/u);      %C�lculos d�o-nos a inclina��o das curvas IS e LM
disp('--------------------------------------------');
       %Reprodu��o e organiza��o do output na janela de comandos
%%%

%% Resolu��o do PVI
subplot(2,2,1);
ode45(@islmdinamico,[t0 tf],[y0 r0]);       %ode45 resolve o PVI e desenha as solu��es de y e de r.
title('Evolu��o do Produto e da taxa de juro');
xlabel('tempo'),ylabel('Produto, taxa de juro');
legend('Produto Y','taxa de juro r');
%%%

%% Defini��o das fun��es de consumo, investimento e importa��es
Cons=C+b*(1-T)*ydot(:,1);
Inv=I-h*ydot(:,2)+jay*ydot(:,1);
Imp=Q+q*ydot(:,1);
%Gr�fico com a evolu��o do produto, consumo, investimento e importa��es
subplot(2,2,2);
plot(t,[ydot(:,1) Cons Inv Imp]);
title('Evolu��o do Produto, Consumo, Investimento e Importa��es');
xlabel('tempo'),ylabel('Produto Y, Cons, Inv, Imp');
legend('Produto Y','Consumo','Investimento','Importa��es');
%%%

%% Diagrama de fase com campos vectoriais 
subplot(2,2,3);
plot(ydot(:,1),ydot(:,2));      %Representa��o gr�fica; Diagrama de fase - Solu��es de y contra r
hold on;
%Campo vectorial - Associa��o de vectores de velocidade a cada ponto das
%solu��es de r e de y associadas � fun��o dynamicISLM utilizando a fun��o
%campovectores.
%Definir valores � grid para o produto yval e para a taxa de juro rval. yval e
%rval t�m que ter a mesma dimens�o pelo que t�m que ser igualmente
%espa�ados. Os passos podem definir-se utilizando a fun�ao linspace.

if min(ydot(:,1))<0     %A margem sobre o m�nimo e m�ximo do produto e taxa de juro serve para melhorar a representa��o gr�fica,
                        %Deve ser definida consoante o sinal do valor
                        %m�nimo do produto e da taxa de juro embora valores
                        %negativos inviabilizem o modelo.
    yval=linspace(min(ydot(:,1))*1.1,max(ydot(:,1))*1.1,20);
else
    yval=linspace(min(ydot(:,1))/1.1,max(ydot(:,1))*1.1,20);
end

if min(ydot(:,2))<0     % "  "
    rval=linspace(min(ydot(:,2))*1.1,max(ydot(:,2))*1.1,20);
else
    rval=linspace(min(ydot(:,2))/1.1,max(ydot(:,2))*1.1,20);
end
campovectorial(@islmdinamico,yval,rval);    %Representa��o gr�fica; Cria��o de campo vectorial(ver fun��o no final).
title('Diagrama de fase com campos vectoriais');
xlabel('Produto');ylabel('taxa de juro');
hold off;
%%%

%% Curva IS e Curva LM - Modelo est�tico no steady state e o ajustamento para
%o steady-state
Y=[-10^10 10^10];     %Definir dois pontos para o produto
R1=(A+a1*Y)/h;                         %Curva IS
R2=(-M/P+L+k*Y)/u;                     %Curva LM
%A assun��o irrealista de um valor negativo para Y trata-se apenas de um
%artif�cio para garantir uma representa��o gr�fica mais "agrad�vel �
%vista".

%Gr�fico
subplot(2,2,4);
plot(ydot(:,1),ydot(:,2));          %Diagrama de fase
axis([min(yval) max(yval) min(rval) max(rval)]);
hold on;
plot(Y,R1,'g',Y,R2,'k');            %IS-LM est�tico
xlabel('Produto'); ylabel('taxa de juro R');
title('Modelo IS-LM - Equil�brio e traject�ria');
legend('traject�ria','IS','LM');
%%%

%% Campos vectoriais, curvas IS e LM, e diagrama de fase para diferentes solu��es e as curvas IS e LM.

%Optar pela representa��o, ou n�o, de um diagrama de fase com n par solu��es de
%y contra as de r dadas n diferentes condi��es iniciais.
ButtonName2 = questdlg('Deseja visualizar diagrama de fase com n par de solu��es y(t) e r(t)? (Nota: s� varia y(0)).', ...
                         'Estabilidade do equil�brio num modelo IS-LM','Sim', 'N�o','N�o');
if strcmp(ButtonName2,'N�o')==1     %Simula��o terminada caso utilizador opte por n�o visualizar o gr�fico.
   return
end
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

%Escolher o n�mero de par solu��es a representar no diagrama de fase
prompt={'N�mero de condi��es iniciais para o produto (para a mesma taxa de juro inicial) '};
numlines=1;
name='Escolher n�mero de par de solu��es y(t) e r(t) a representar';
defaultanswer={'1'};     
answer4=inputdlg(prompt,name,numlines,defaultanswer,options); 
if isempty(answer4)     
    return                  %Simula��o concluida se utilizador decidir "cancelar"
end

nrtraj=str2double(answer4{1});       %nrtraj - N�mero de condi��es iniciais y(0)  

figure;                 
hold on;

y01=linspace(1,y0,nrtraj);      %r01 linearmente espa�ado por nrtraj(nr de condi��es iniciais) pontos, torna gr�fico mais agrad�vel � vista.
AA=zeros(1,nrtraj); BB=AA;
CC=AA; DD=AA;

for i= 1:nrtraj
    [~,ydotv]=ode45(@islmdinamico,[t0 tf],[y01(i) r0]);
    AA(i)=max(ydotv(:,1)); BB(i)=min(ydotv(:,1));
    CC(i)=max(ydotv(:,2)); DD(i)=min(ydotv(:,2));
end                 %As matrizes t�m m�ximos e m�nimos da taxa de juro e do produto de todas as solu��es de y(t) e r(t) tomadas no seu conjunto, para definir as margens do campo vectorial.
y2val=linspace(min(BB),max(AA),20);
r2val=linspace(min(DD),max(CC),20);
campovectorial(@islmdinamico,y2val,r2val);      %Campos vectoriais
plot(Y, R1,'g',Y, R2,'k')           %Curvas IS  e LM
for y01=linspace(1,y0,nrtraj)       %Aqui repete-se o comando em vez de ter feito os plot's todos anteriormente para n�o afectar a legendagem do gr�fico.
    [~,ydotv]=ode45(@islmdinamico,[t0 tf],[y01 r0]);
    plot(ydotv(:,1),ydotv(:,2));        %V�rias traject�rias no diagrama de fase.          
end
title('IS-LM: Diagrama de fase para diferentes valores iniciais da taxa de juro nominal');
xlabel('Produto');ylabel('Taxa de juro nominal')
legend('Vectores - Direc��es','IS','LM','Traject�rias');
hold off;
%%%

%% Fun��es
%Sistema de duas equa��es diferenciais representativo do modelo IS-LM
function [ ydot ] = islmdinamico( ~,y )   
global alpha C I G X Q b T jay g q h beta M P L k u
%Defini��o do modelo IS-LM din�mico (pre�os fixos, ex�genos)
a1=b*(1-T)+jay+g-q-1;       %Inverso do sim�trico do multiplicador keynesiano a1
A=C+G+I+X-Q;                %Despesa aut�noma
ydot=zeros(2,1);        
ydot(1)=alpha*(A+a1*y(1)-h*y(2));       %Din�mica do mercado de bens e servi�os
ydot(2)=beta*(-M/P+L+k*y(1)-u*y(2));    %Din�mica do mercado monet�rio
%Em que y(1) e y(2) s�o, respectivamente, as solu��es do produto real e da taxa
%de juro nominal.

%% Campo vectorial - Associa��o de vectores de velocidade u,v a cada ponto
% x,y e desenho do gr�fico com os vectores
function campovectorial(func,y1val,r1val)

n1=length(y1val);
n2=length(r1val);
%yp1 e yp2 s�o os vectores associados a yval e rval, respectivamente.
yp1=zeros(n2,n1);       %Indiferente pois y1val e y2val t�m que ser igualmente espa�ados
yp2=zeros(n2,n1);
for i=1:n1
  for j=1:n2
    ypv = feval(func,0,[y1val(i);r1val(j)]);        %yp1 e yp2 s�o os vectores associados a yval e rval, respectivamente.
    yp1(j,i) = ypv(1);
    yp2(j,i) = ypv(2);
  end
end
quiver(y1val,r1val,yp1,yp2,'c','LineWidth',1.4);        %Representa��o gr�fica; Desenha setas nos vectores de velocidade
axis tight;