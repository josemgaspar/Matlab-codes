%% Modelo IS-LM dinâmico
%Função que realiza uma análise ao modelo IS-LM (economia fechada ou grande
%economia aberta) dinâmico. Incorpora a utilização de um solver ode45 para
%sistemas de equações diferenciais para representar as soluções do produto
%e da taxa de juro ao longo do tempo (Problema de valor inicial); representa as soluções de y(t)
%contra as de r(t) em diagramas de fase. Para os diagramas de fase, cria
%ainda um campo vectorial com vectores de velocidade que permite inferir e
%tornar intuitivas as direcções das trajectórias nos diagramas de fase. Representa ainda
%a evolução das principais variáveis agregadas macroeconómicas e faz uma
%análise à estabilidade dos pontos de quilíbrio do sistema e à inclinação das curvas IS e LM.
%O output gerado são valores concretos para algumas das variáveis
%mencionadas e 5 gráficos, sendo os primeiros 4 agrupados numa figura e o
%outro numa 2ª.

function modeloISLM
clc; clear all; close all;      %Fechar todas as figuras abertas e limpar a janela de comandos
global alpha C I G X Q b T jay g q h beta M P L k u
disp('------------------------');
disp('Dinâmica do modelo IS-LM'),
disp('------------------------');
%%Definição de variáveis
%Flexibilizar a box gerada pelo comando inputdlg
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

%% A atribuição de valores para as variáveis é feita através de uma box de diálogo.

%Variáveis que não afectam estabilidade do equilíbrio:
prompt = {'Componente autónoma do consumo C=','Componente autónoma do Investimento I=',...
    'Componente autónoma dos gastos públicos G=','Exportações X=','Componente autónoma das importações Q=',...
    'Oferta nominal de moeda M=','Nível geral de preços P=','Procura autónoma de moeda L='};
name = 'Parâmetros do modelo IS-LM (não afectam estabilidade do equilíbrio)';
numlines = 1;
defaultanswer = {'6','25','15','6','9','30','1','3'};     %Valores atribuídos por defeito - Eqº estável.
answer1 = inputdlg(prompt,name,numlines,defaultanswer,options);       %Faz aparecer a box para inserir os valores.

if isempty(answer1)     %Sai do programa se não forem atribuídos valores (escolher cancelar).
    clc;
    return
end

C=str2double(answer1{1});       %Consumo autónomo
I=str2double(answer1{2});       %Investimento autónomo
G=str2double(answer1{3});       %Despesa Pública autónoma
X=str2double(answer1{4});       %Exportações
Q=str2double(answer1{5});       %Importação autónoma
M=str2double(answer1{6});       %Oferta nominal de moeda
P=str2double(answer1{7});       %Nível geral de preços
L=str2double(answer1{8});       %Procura autónoma por moeda

%Variáveis que afectam a estabilidade do sistema
prompt={'Parâmetro da equação diferencial dy/dt=alpha(A+a1*y(t)-hr(t)  alpha=',...
    'Parâmetro da equação diferencial dr/dt=beta(-M/P+L+ky(t)-ur(t))  beta=','Taxa marginal de imposto T=',...
    'Propensão marginal a consumir b=','Sensibilidade do investimento ao produto real jay=',...
    'Sensibilidade do investimento à taxa de juro nominal h=','Sensibilidade da procura real de moeda face ao produto real k=',...
    'Sensibilidade da procura real de moeda face à taxa de juro nominal u=','Propensão marginal à despesa pública g=',...
    'Propensão marginal à importação q='};
name='Parâmetros do modelo IS-LM (afectam a estabilidade do equilíbrio)';
numlines=1;
defaultanswer={'0.5','0.4','0.25','0.7','0.3','1.3','0.25','0.3','0.05','0.1'};        %Valores atribuídos por defeito - Eqº estável.
answer2=inputdlg(prompt,name,numlines,defaultanswer,options);       % ""

if isempty(answer2)     %""
    clc;
    return
end

alpha=str2double(answer2{1});       %Parâmetro alpha da equação diferencial dy/dt
beta=str2double(answer2{2});        %Parâmetro beta da equação diferencial dr/dt
T=str2double(answer2{3});           %Taxa de imposto marginal
b=str2double(answer2{4});           %Propensão marginal ao consumo
jay=str2double(answer2{5});         %Sensibilidade do investimento ao produto real
h=str2double(answer2{6});           %Sensibilidade do investimento à taxa de juro nominal
k=str2double(answer2{7});           %Sensibilidade da procura real de moeda face ao produto real
u=str2double(answer2{8});           %Sensibilidade da procura real de moeda face à taxa de juro nominal
g=str2double(answer2{9});           %Propensão marginal à despesa pública
q=str2double(answer2{10});          %Propensão marginal à importação

%Valor inicial e final para o tempo e condições iniciais y(0) e r(0)
prompt={'Valor inicial do produto real  y(0)=','Valor inicial da taxa de juro nominal  r(0)=',...
    'tempo inicial  t0=','tempo final  tf='};
name='PVI: Condições iniciais';
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

a1=b*(1-T)+jay+g-q-1;       %Inverso do simétrico do multiplicador keynesiano a1
A=C+G+I+X-Q;                %Despesa autónoma
%%%

%% Apresentação do modelo e parâmetros usados
 %Mostrar na janela de comandos os valores escolhidos para os parâmetros
fprintf('\nModelo dinâmico IS-LM:\n                           dy/dt=alpha(A+a1*y(t)-hr(t)\n                           dr/dt=beta(-M/P+L+ky(t)-ur(t))\n');
fprintf('\nCurvas IS e LM:\n                    (IS)r=(A+a1*y)/h\n                    (LM)r=(k*y-M/P+L)/u\n ');
fprintf('\ny: Produto real\nr: Taxa de juro nominal\nA: Despesa autónoma = C+G+I+X-Q = %g\na1: b*(1-T)+jay+g-q-1 = %g\n',...
    A,a1);
fprintf('\nParâmetros usados:\nC = %g; I =  %g; G = %g; X = %g; Q = %g; M = %g; P = %g; L = %g; alpha = %g; beta = %g;\n',...
    C,I,G,X,Q,M,P,L,alpha,beta);
fprintf('T = %g; b = %g; jay = %g; h = %g; k = %g; u = %g; g = %g; q = %g\n------------------------------------------------------------\n',...
    T,b,jay,h,k,u,g,q);                                                    

%O programa dá erro se o produto de equilíbrio for negativo
if k*h-a1*u<0
    errordlg('Produto de equilíbrio é negativo','Erro da Função');
end

%% PVI: Problema de valor inicial
[t,ydot]=ode45(@islmdinamico,[t0 tf],[y0 r0]);      %O solver ode45 atribuir a ydot os 
%valores das soluções das equações diferenciais definidas em
%@islmdinamico(ver no final).

%Produto e/ou taxa de juro não devem ser negativos em nenhum momento do
%tempo
if all(ydot(:) > 0)        %Se todos os valores do produto e da taxa de juro forem positivas, continua.
else
    %Cria uma box em que é dado a escolher continuar com a função ou não,
    %caso y ou r negativos em algum tempo
    ButtonName = questdlg('O valor escolhido para os parâmetros implica níveis do produto e/ou taxa de juro negativos em algum momento do tempo. Continuar simulação?', ...
                         'Atenção!','Sim', 'Não','Não');
                     fprintf('\n!!!Produto e/ou taxa de juro não podem ser inferiores ou iguais a zero em nenhum momento do tempo:\nSimulação meramente ilustrativa da dinâmica do modelo. \n\n')
                     if strcmp(ButtonName,'Não')==1     %Cancela o programa
                         fprintf('\nSimulação cancelada: Produto e/ou taxa de juro não podem ser inferiores ou iguais a zero em nenhum momento do tempo\n\n');
                             return
                    end
end
%%%


%% Cálculo do ponto de equilíbrio do sistema de equações diferenciais(steady-state).
options=optimset('Display','off');      %Suprimir informações do solver fsolve
[pe]=fsolve(@(y) islmdinamico(t,y),[y0;r0],options);        %Calcula y e r de equilíbrio (solução local) dado dy/dt= e dr/dt=0
%e utilizando y0 e r0 como condições iniciais.
disp(' ');
disp('Steady-State:');
fprintf('\nProduto de equilíbrio    Taxa de juro de equilíbrio\n      %f                 %f\n',pe(1,1),pe(2,1));        %Retorna o ponto de eqº

%
%% Estudo da estabilidade do ponto de equilíbrio
% Estudo da estabilidade em função dos valores próprios da matriz jacobiana
Jac=[alpha.*a1 -alpha*h;beta*k -beta*u];        %Matriz jacobiana no ponto de equilíbrio.
disp(' ');
disp('Valores próprios da matriz jacobiana');
disp(eig(Jac));        %Mostra os val pp da matriz jacobiana


if real(eig(Jac))<0 
    if isreal(eig(Jac))==0  %Para ser estável, ambos os valores próprios têm que ter parte real negativa
    disp('O equilíbrio é um foco estável.');
    else
        disp('O equilíbrio é um nó estável.')
    end
elseif real(eig(Jac))>0 
    if isreal(eig(Jac))==0
    disp('O equilíbrio é um foco instável.')
    elseif isreal(eig(Jac))==1
    disp('O equilíbrio é um nó instável.')
    end
elseif real(eig(Jac))==0 
    if isreal(eig(Jac))==0
    disp('O equilíbrio é um ponto centro estável. Existem soluções periódicas');
    end
else    
    disp('O equilíbrio é instável.')        %Não é possivel haver pontos sela, pois det(Jac)>0, já que queremos um y de equilíbrio positivo.
end

if a1<0
    fprintf('\nEstabilidade do equilíbrio devido a:\nMultiplicador keynesiano maior que zero (a1<0)\n');
elseif a1>0 && alpha*a1<beta*u      %Relação entre a1 e Multiplicador e consequências dos seus sinais (e/ou do traço da jacobiana) sobre a estabilidade.
    fprintf('\nEstabilidade do equilíbrio devido a:\nMultiplicador keynesiano menor que zero (a1>0) e tr(Jac) < 0\n(nota:a1=b*(1-T)+jay+g-q-1\n');
else 
    fprintf('\nInstabilidade do equilíbrio devido a:\nMultiplicador keynesiano menor que zero (a1>0) e tr(Jac) > 0\n(nota:a1=b*(1-T)+jay+g-q-1\n');
end

fprintf('\nValor de outras variáveis agregadas no steady-state:\n');     %Organização e reprodução d output na janela de comando

%Consumo, Investimento, Despesa pública, Exportações e Importações no
%steady-state: Output na janela de comandos
Ceq=C+b*(1-T)*pe(1,1);
Ieq=I-h*pe(2,1)+jay*pe(1,1);
Imp=Q+q*pe(1,1);
Geq=G+g*pe(1,1);
fprintf('\nConsumo    Investimento    Despesa Pública     Exportações      Importações\n %.3f        %.3f          %.3f             %.3f            %.3f\n',...
    Ceq,Ieq,Geq,X,Imp);
disp('--------------------------------------------------------------------------------');
%%%

%% Multiplicadores do modelo IS-LM estático
mult1=-1/a1;     %Multiplicador keynesiano
mult2=-1/(a1-h*k/u);        %Multiplicador dos gastos públicos
mult3=-(k/u)/(a1-h*k/u);      %Multiplicador monetário
mult4=b*pe(1,1)/(a1-h*k/u);      %Multiplicador da taxa de imposto
fprintf('\nMultiplicadores (Modelo IS-LM estático):\n \nMultiplicador keynesiano (-1/a1):  %f\nMultiplicador dos gastos públicos: %f \nMultiplicador monetário:           %f \nMultiplicador da taxa de imposto: %f\n',...
    mult1,mult2,mult3,mult4);
disp('--------------------------------------------');
disp(' ');
disp('Inclinação das curvas IS e LM:');     
disp(' ');
if a1<0         %A inclinação pode ser definida em função do multiplicador keynesiano (relação com a1)
    disp('(Multiplicador keynesiano maior que zero (a1<0): IS com inclinação negativa)');
elseif a1==0
    disp('(Multiplicador keynesiano indeterminado (a1=0): IS com inclinação nula; IS horizontal)');
else
    disp('(Multiplicador keynesiano menor que zero (a1>0): IS com inclinação positiva)');
end
fprintf('\nInclinação IS (a1/h)    Inclinação LM (k/u)\n      %f              %f\n',a1/h,k/u);      %Cálculos dão-nos a inclinação das curvas IS e LM
disp('--------------------------------------------');
       %Reprodução e organização do output na janela de comandos
%%%

%% Resolução do PVI
subplot(2,2,1);
ode45(@islmdinamico,[t0 tf],[y0 r0]);       %ode45 resolve o PVI e desenha as soluções de y e de r.
title('Evolução do Produto e da taxa de juro');
xlabel('tempo'),ylabel('Produto, taxa de juro');
legend('Produto Y','taxa de juro r');
%%%

%% Definição das funções de consumo, investimento e importações
Cons=C+b*(1-T)*ydot(:,1);
Inv=I-h*ydot(:,2)+jay*ydot(:,1);
Imp=Q+q*ydot(:,1);
%Gráfico com a evolução do produto, consumo, investimento e importações
subplot(2,2,2);
plot(t,[ydot(:,1) Cons Inv Imp]);
title('Evolução do Produto, Consumo, Investimento e Importações');
xlabel('tempo'),ylabel('Produto Y, Cons, Inv, Imp');
legend('Produto Y','Consumo','Investimento','Importações');
%%%

%% Diagrama de fase com campos vectoriais 
subplot(2,2,3);
plot(ydot(:,1),ydot(:,2));      %Representação gráfica; Diagrama de fase - Soluções de y contra r
hold on;
%Campo vectorial - Associação de vectores de velocidade a cada ponto das
%soluções de r e de y associadas à função dynamicISLM utilizando a função
%campovectores.
%Definir valores à grid para o produto yval e para a taxa de juro rval. yval e
%rval têm que ter a mesma dimensão pelo que têm que ser igualmente
%espaçados. Os passos podem definir-se utilizando a funçao linspace.

if min(ydot(:,1))<0     %A margem sobre o mínimo e máximo do produto e taxa de juro serve para melhorar a representação gráfica,
                        %Deve ser definida consoante o sinal do valor
                        %mínimo do produto e da taxa de juro embora valores
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
campovectorial(@islmdinamico,yval,rval);    %Representação gráfica; Criação de campo vectorial(ver função no final).
title('Diagrama de fase com campos vectoriais');
xlabel('Produto');ylabel('taxa de juro');
hold off;
%%%

%% Curva IS e Curva LM - Modelo estático no steady state e o ajustamento para
%o steady-state
Y=[-10^10 10^10];     %Definir dois pontos para o produto
R1=(A+a1*Y)/h;                         %Curva IS
R2=(-M/P+L+k*Y)/u;                     %Curva LM
%A assunção irrealista de um valor negativo para Y trata-se apenas de um
%artifício para garantir uma representação gráfica mais "agradável à
%vista".

%Gráfico
subplot(2,2,4);
plot(ydot(:,1),ydot(:,2));          %Diagrama de fase
axis([min(yval) max(yval) min(rval) max(rval)]);
hold on;
plot(Y,R1,'g',Y,R2,'k');            %IS-LM estático
xlabel('Produto'); ylabel('taxa de juro R');
title('Modelo IS-LM - Equilíbrio e trajectória');
legend('trajectória','IS','LM');
%%%

%% Campos vectoriais, curvas IS e LM, e diagrama de fase para diferentes soluções e as curvas IS e LM.

%Optar pela representação, ou não, de um diagrama de fase com n par soluções de
%y contra as de r dadas n diferentes condições iniciais.
ButtonName2 = questdlg('Deseja visualizar diagrama de fase com n par de soluções y(t) e r(t)? (Nota: só varia y(0)).', ...
                         'Estabilidade do equilíbrio num modelo IS-LM','Sim', 'Não','Não');
if strcmp(ButtonName2,'Não')==1     %Simulação terminada caso utilizador opte por não visualizar o gráfico.
   return
end
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

%Escolher o número de par soluções a representar no diagrama de fase
prompt={'Número de condições iniciais para o produto (para a mesma taxa de juro inicial) '};
numlines=1;
name='Escolher número de par de soluções y(t) e r(t) a representar';
defaultanswer={'1'};     
answer4=inputdlg(prompt,name,numlines,defaultanswer,options); 
if isempty(answer4)     
    return                  %Simulação concluida se utilizador decidir "cancelar"
end

nrtraj=str2double(answer4{1});       %nrtraj - Número de condições iniciais y(0)  

figure;                 
hold on;

y01=linspace(1,y0,nrtraj);      %r01 linearmente espaçado por nrtraj(nr de condições iniciais) pontos, torna gráfico mais agradável à vista.
AA=zeros(1,nrtraj); BB=AA;
CC=AA; DD=AA;

for i= 1:nrtraj
    [~,ydotv]=ode45(@islmdinamico,[t0 tf],[y01(i) r0]);
    AA(i)=max(ydotv(:,1)); BB(i)=min(ydotv(:,1));
    CC(i)=max(ydotv(:,2)); DD(i)=min(ydotv(:,2));
end                 %As matrizes têm máximos e mínimos da taxa de juro e do produto de todas as soluções de y(t) e r(t) tomadas no seu conjunto, para definir as margens do campo vectorial.
y2val=linspace(min(BB),max(AA),20);
r2val=linspace(min(DD),max(CC),20);
campovectorial(@islmdinamico,y2val,r2val);      %Campos vectoriais
plot(Y, R1,'g',Y, R2,'k')           %Curvas IS  e LM
for y01=linspace(1,y0,nrtraj)       %Aqui repete-se o comando em vez de ter feito os plot's todos anteriormente para não afectar a legendagem do gráfico.
    [~,ydotv]=ode45(@islmdinamico,[t0 tf],[y01 r0]);
    plot(ydotv(:,1),ydotv(:,2));        %Várias trajectórias no diagrama de fase.          
end
title('IS-LM: Diagrama de fase para diferentes valores iniciais da taxa de juro nominal');
xlabel('Produto');ylabel('Taxa de juro nominal')
legend('Vectores - Direcções','IS','LM','Trajectórias');
hold off;
%%%

%% Funções
%Sistema de duas equações diferenciais representativo do modelo IS-LM
function [ ydot ] = islmdinamico( ~,y )   
global alpha C I G X Q b T jay g q h beta M P L k u
%Definição do modelo IS-LM dinâmico (preços fixos, exógenos)
a1=b*(1-T)+jay+g-q-1;       %Inverso do simétrico do multiplicador keynesiano a1
A=C+G+I+X-Q;                %Despesa autónoma
ydot=zeros(2,1);        
ydot(1)=alpha*(A+a1*y(1)-h*y(2));       %Dinâmica do mercado de bens e serviços
ydot(2)=beta*(-M/P+L+k*y(1)-u*y(2));    %Dinâmica do mercado monetário
%Em que y(1) e y(2) são, respectivamente, as soluções do produto real e da taxa
%de juro nominal.

%% Campo vectorial - Associação de vectores de velocidade u,v a cada ponto
% x,y e desenho do gráfico com os vectores
function campovectorial(func,y1val,r1val)

n1=length(y1val);
n2=length(r1val);
%yp1 e yp2 são os vectores associados a yval e rval, respectivamente.
yp1=zeros(n2,n1);       %Indiferente pois y1val e y2val têm que ser igualmente espaçados
yp2=zeros(n2,n1);
for i=1:n1
  for j=1:n2
    ypv = feval(func,0,[y1val(i);r1val(j)]);        %yp1 e yp2 são os vectores associados a yval e rval, respectivamente.
    yp1(j,i) = ypv(1);
    yp2(j,i) = ypv(2);
  end
end
quiver(y1val,r1val,yp1,yp2,'c','LineWidth',1.4);        %Representação gráfica; Desenha setas nos vectores de velocidade
axis tight;