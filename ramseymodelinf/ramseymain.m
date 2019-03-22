%Problema de optimização dinâmica - Crescimento óptimo
%Modelo de Ramsey-Cass-Koopmans com progresso técnico
%em equilíbrio centralizado com
%horizonte infinito.

function ramseymain

clear all
close all

global rho delta alpha sigma gx n guessbvp NMax RelTol AbsTol kss k0 c0...
    p A0 N0 Q0

%---------------------------------------------------------------%
%%Parâmetros base
%---------------------------------------------------------------%
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';

prompt = {'População atual, Q_0=',...
    'Taxa de crescimento do progresso técnico A, gx =',...
    'Nível tecnológico, A_0 =',...
    'Horas de trabalho atuais, N_0=',...
    'Taxa de crescimento da oferta de trabalho n =',...
    'Taxa de depreciação do capital físico \delta =',...
    'Elasticidade do produto ao capital \alpha =',...
    'Taxa de desconto (factor de preferência pelo tempo) \rho =',...
    'Elasticidade intertemporal do consumo \sigma =',...
    'k_0 em fracção do steady-state, k(0) / kss =','t(0) =','t(T) ='};
name = 'Parâmetros do modelo';
numlines = 1;
defaultanswer = {'1','0.02','1','1','0.04','0.2','1/3','0.08','3','0.5','0','100'};     
answer1 = inputdlg(prompt,name,numlines,defaultanswer,options);     

if isempty(answer1)    
    clc;
    return
end

Q0 = eval(answer1{1});
gx = eval(answer1{2});
A0 = eval(answer1{3});
N0 = eval(answer1{4});
n = eval(answer1{5});       %
delta = eval(answer1{6});       %
alpha = eval(answer1{7});       
rho = eval(answer1{8});     
sigma = eval(answer1{9});      
p = eval(answer1{10});       
t0 = eval(answer1{11});
tf = eval(answer1{12});

kss = (alpha/(delta+rho+gx/sigma))^(1/(1-alpha));       % k de steady-state
css = kss^alpha-(delta+n+gx)*kss;                       % c de steady-state                                         % t(0) e t(T)
k0 = kss*p;                                        % k(0) = k0 > 0
c0 = k0^alpha-(n+delta+gx)*k0;                          % c(0)= c0
%---------------------------------------------------------------%


%-------------------------------------------------------------------%
%Escolher o método para resolver o problema de optimização dinâmica
%-------------------------------------------------------------------%
ButtonName2 = questdlg...
('Modelo de Ramsey-Cass-Koopmans - Método Directo... (GPOPS) vs Indirecto (BVP5C)', ...
                         'Modelo de Ramsey','GPOPS', 'BVP5C','Cancelar','Cancelar');

                     
%Resolução directa - GPOPS                     
if strcmp(ButtonName2,'GPOPS')==1     

x = cputime;
ramseyGPOPS;
w = cputime-x;
t = solution.time;
k = solution.state;
c = solution.control;
lambda = -solution.costate; %Queremos um multiplicador lagrangeano com
                            %valores positivos

%Resolução indirecta - BVP5C
elseif strcmp(ButtonName2,'BVP5C')== 1 

x = cputime;
ramseyBVP5C;
w = cputime-x;
t = linspace(t0,tf);
y = deval(sol,t);
k = y(1,:)';
c = y(2,:)';
t = t';
       
else     
    return  % Cancelar programa   
end
%-------------------------------------------------------------------%


%---------------------------------------------------------------%        
%Desenhar figuras
%---------------------------------------------------------------%

figure(1);
pp = plot(t,k,'-o');
set(pp,'LineWidth',1.5);
set(gca,'FontName','Times','FontSize',16);
xl = xlabel('$t$','Interpreter','latex');
yl = ylabel('$k_a(t)$','Interpreter','latex');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
title('Capital por unidade eficiente de trabalho');
grid on;

figure(2);
pp = plot(t,c,'-o');
set(pp,'LineWidth',1.5);
set(gca,'FontName','Times','FontSize',16);
xl = xlabel('$t$','Interpreter','latex');
yl = ylabel('$c_a(t)$','Interpreter','latex');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
title('Consumo por unidade eficiente de trabalho');
grid on;

if exist('lambda','var') == 0
	lambda = costate(c,rho,sigma,n,gx,t);
else
      
end

figure(3);
pp = plot(t,lambda,'-o');
set(pp,'LineWidth',1.5);
set(gca,'FontName','Times','FontSize',16);
xl = xlabel('$t$','Interpreter','latex');
yl = ylabel('$\lambda(t)$','Interpreter','latex');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
title('Multiplicador Lagrange');
grid on; 

%Figura 4: Definir variáveis para o diagrama de fase

ButtonName8 = questdlg...
('Visualizar diagrama de fase?', ...
                         'Diagrama de fase','Sim', 'Não','Não');
                     
if strcmp(ButtonName8,'Sim')== 1
kval = linspace(min(k)*0.8,max(k)*1.2,18);  
cval = linspace(min(c)*0.8,max(c)*1.2,18);       
kp = linspace(min(k)*0.8,max(k)*1.2,150);
c_dk_0 = kp.^alpha-(delta+n+gx).*kp; % dk/dt=0 => c(t)=(...)

figure(4);
pp = plot(k,c);
set(pp,'LineWidth',1.5);
set(gca,'FontName','Times','FontSize',16);
xl = xlabel('$k(t)$','Interpreter','latex');
yl = ylabel('$c(t)$','Interpreter','latex');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
grid on;
hold on;
ppp = plot(kp,c_dk_0,'k');
set(ppp,'LineWidth',1.5)      

if min(c_dk_0) < min(c)
        line([kss kss],[min(c_dk_0)*0.8 max(c)*1.2],[1 1],...
           'Color','g','LineWidth',1.5) % Linha dc/dt = 0
else
        line([kss kss],[min(c)*0.8 max(c)]*1.2,[1 1],...
           'Color','g','LineWidth',1.5)
end  

ypv = [gradient(k) gradient(c)]; %tirar, se função não for suave ou poucos pontos
u = ypv(:,1);
v = ypv(:,2);
%Desenha setas nos vectores de velocidade
quiver(k,c,u,v,'b','LineWidth',1.5);  
axis tight;
campovectorial(@eqdinamicas,kval,cval);
leg = legend('$c_{(k)}$','$\frac{dc_a}{dt}=0$','$\frac{dk_a}{dt}=0$');
set(leg,'Interpreter','latex');
title('Diagrama de Fase');
else
end
%---------------------------------------------------------------%


%---------------------------------------------------------------%
%Variaveis em nível, per capita, e por hora de trabalho
%---------------------------------------------------------------%
q = n; %Taxa de crescimento da população determina a tx. crescimento
%oferta de trabalho
Q = Q0.*exp(q.*t); %Evolução da população
ya = gdp(k);    %produto por unidade de trabalho eficiente
s = savings(ya,c); %poupança por unidade eficiente de trabalho
yph = ya.*A0.*exp(gx.*t); %produto por hora de trabalho
cph = c.*A0.*exp(gx.*t);    %consumo por hora de trabalho
kph = k.*A0.*exp(gx.*t);    %capital por hora de trabalho
sph = s.*A0.*exp(gx.*t);    %poupança por hora de trabalho
Y = yph.*N0.*exp(n.*t);     %Produto real
K = kph.*N0.*exp(n.*t);     %Stock de capital
C = cph.*N0.*exp(n.*t);     %Consumo
S = sph.*N0.*exp(n.*t);     %Poupança
ypc = Y./Q;                 %Pib per capita
cpc = C./Q;                 %Consumo per capita
kpc = K./Q;                 %Capital per capita
spc = S./Q;                 %Poupança per capita

ButtonName3 = questdlg...
('Gráfico das variáveis per capita?', ...
                         'Gráficos','Sim', 'Não','Não');
                     
if strcmp(ButtonName3,'Sim')==1

figure
pp = plot(t,ypc,'y',t, spc,'r',t, kpc, 'c', t, cpc, 'b');
set(gca,'FontName','Times','FontSize',16);
yl = ylabel('Variáveis per capita');
xl = xlabel('t','Interpreter','latex');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
title('Variáveis per capita');
legend('Produto real pc','Poupança pc','Capital pc', 'Consumo pc');
set(pp,'LineWidth',1.5);
else
end
%---------------------------------------------------------------%

%---------------------------------------------------------------%
%BALANCED GROWTH PATH
%---------------------------------------------------------------% 

%Em termos continuos
BGPka = kgr(@eqdinamicas,k,c);
dya = (0.5.*k.^0.5).*kdot; %dy/dt = (dy/dk)*(dk/dt)
BGPya = dya./ya;
BGPca = cgr(@eqdinamicas,k,c);
BGPkph = BGPka+gx;
BGPyph = BGPya+gx;
BGPcph = BGPca+gx;
BGPkpc = BGPkph;
BGPcpc = BGPcph;
BGPypc = BGPyph;
BGPK = BGPkph+n;
BGPY = BGPyph+n;
BGPC = BGPypc+n;


ButtonName4 = questdlg...
('Visualizar gráfico das taxas de crescimento?', ...
                         'Gráfico','Sim', 'Não','Não');
                     
if strcmp(ButtonName4,'Sim')== 1

figure
pp = plot(t,BGPypc,'b',t,BGPcpc,'r',t,BGPkpc,'k');
set(gca,'FontName','Times','FontSize',16);
set(pp,'LineWidth',1.5);
labels = get(gca,'yticklabel');
labels_modif=[num2str(100*str2num(labels)) ones(length(labels),1)*'%'];  %#ok<ST2NM>
set(gca,'yticklabel',labels_modif);
title('Balanced Growth Path','Interpreter','latex');
yl = ylabel('Taxas de crescimento','Interpreter','latex');
xl = xlabel('t','Interpreter','latex');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
legend('Tx. cresc. Produto pc','Tx. cresc. Consumo pc','Tx. cresc. Capital pc');
else
end
%---------------------------------------------------------------%

%---------------------------------------------------------------%
%Output na janela de comandos:
%---------------------------------------------------------------%
clc;

disp('----------------------------------------------------');
disp('Modelo de Ramsey-Kass-Coopmans com progresso técnico'),
disp('----------------------------------------------------');
fprintf('\nParâmetros usados:\nrho = %g; delta =  %g; alpha = %g; \n',...
    rho,delta,alpha);
fprintf('sigma = %g; gx = %g; n = %g; k0 = %g; c0 = %g;\n',...
     sigma,gx,n,k0,c0);
fprintf('A0 = %g; N0 = %g; Q0 = %g\n',A0,N0,Q0);     
fprintf('\n------------------------------------------------------------\n')

if BGPka(l) <= 1e-3
    disp('Atingiu o steady-state.')
    disp(['kss: ',num2str(kss),'  css: ',num2str(css)]);
else
    disp('Não atingiu o steady-state!')
end
disp('');
disp('Balanced growth path:')
disp(['Taxa de crescimento ya: ',num2str(BGPya(l)*100),'%']);
disp(['Taxa de crescimento produto per capita: ',...
    num2str(BGPypc(l)*100),'%']);
disp(['Taxa de crescimento produto real: ',num2str(BGPY(l)*100),'%']);
disp(['Cpu time: ',num2str(w)]);
%---------------------------------------------------------------%

%---------------------------------------------------------------%
%Soluçao analitica?
%---------------------------------------------------------------%
ButtonName5 = questdlg...
('Incluir solução analítica ou aproximação?', ...
                         'Sol. analítica','Sim', 'Não','Não');   
                       
if strcmp(ButtonName5,'Sim') == 1   

prompt = {'Solução capital (ksol=)','Solução consumo (csol=...)'};
name = 'Expressão de solução analítica';
numlines = 1;
%Nota: Esta solução por defeito só funciona com sigma = 1/alpha
defaultanswer =...
    {['((kss.^(1-alpha)+(k0.^(1-alpha)-kss.^(1-alpha)).*exp(-(1-alpha).*'...
    '((rho+delta+gx/sigma)/alpha).*t)).^(1/(1-alpha)))'],...
    ['(((rho-n-gx*(1-1/sigma))+(n+delta+gx)*(1-alpha))'...
    '.*(((kss.^(1-alpha)+(k0.^(1-alpha)-kss.^(1-alpha)).*exp(-(1-alpha)'...
    '.*((rho+delta)/alpha).*t)).^(1/(1-alpha)))))/alpha;']};     %Valores atribuídos por defeito - Eqº estável.
answer2 = inputdlg(prompt,name,numlines,defaultanswer,options);      


if isempty(answer2)
else


    
if strcmp(answer2,defaultanswer') == 1 & sigma ~= 1/alpha %#ok<AND2>


    f = warndlg('Esta solução só é válida para sigma = 1/alpha!','Atenção!');
    waitfor(f);


else
ksol = eval(answer2{1});
csol = eval(answer2{2});
lambdasol = costate(csol,rho,sigma,n,gx,t);
        
figure(1);
hold on;
plot(t,ksol,'r','LineWidth',1.5);
legend('Aproximação','sol. analítica');
hold off;

figure(2);
hold on;
plot(t,csol,'r','LineWidth',1.5);
legend('Aproximação','sol. analítica');
hold off;
figure(3)
hold on;
plot(t,lambdasol,'r','LineWidth',1.5);
legend('Aproximação','sol. analítica');
hold off;
assignin('base','lambdaanalitico',lambdasol);
assignin('base','kanalitico',ksol);
assignin('base','canalitico',csol);
end
end
%---------------------------------------------------------------%

%---------------------------------------------------------------%
% Adicionar as variáveis ao "Workspace"
%---------------------------------------------------------------%
assignin('base','t',t);     %vetor do tempo
assignin('base','ca',c);        %C por unidade eficiente de trabalho
assignin('base','ka',k);        %K por unidade eficiente de trabalho
assignin('base','kss',kss);     %ka de steady-state
assignin('base','css',css);     %ca de steady-state
assignin('base','ya',ya);       %Y po unidade eficiente de trabalho
assignin('base','Capitalpc',kpc);   %K per capita
assignin('base','Consumopc',cpc');  %C per capita
assignin('base','Produtopc',ypc);   %Y per capita
assignin('base','sa',s);        %poupança por unidade eficiente de trabalho
assignin('base','Poupancapc',spc);  %S per capita
assignin('base','Produtoph',yph);   %Y por hora de trabalho
assignin('base','Capitalph',kph);   %K por hora de trabalho
assignin('base','Consumoph',cph);   %C por hora de trabalho
assignin('base','Poupancaph',sph);  %C por hora de trabalho
assignin('base','Produto',Y);       %Produto
assignin('base','Capital',K);       %Capital
assignin('base','Poupanca',S);      %Poupança
assignin('base','Consumo',C);       %Consumo
assignin('base','BGPya',BGPya);     %tx crescimento ya
assignin('base','BGPypc',BGPypc);   %tx crescimento ypc
assignin('base','BGPY',BGPY);       %       "       Y
assignin('base','BGPkpc',BGPkpc);   %       "       kpc
assignin('base','BGPcpc',BGPcpc);   %       "       cpc
assignin('base','BGPK',BGPK);       %       "       K
assignin('base','BGPC',BGPC);       %       "       C
assignin('base','Coestado',lambda); %Multiplicador dinamico, ou co-estado
end
%---------------------------------------------------------------%

%---------------------------------------------------------------%
%Comparar os métodos. Só funciona se primeiro foi escolhido GPOPS
%---------------------------------------------------------------%
if strcmp(ButtonName2,'GPOPS')==1   
    
ButtonName6 = questdlg...
('Comparar GPOPS com BVP5C?', ...
                         'Comparar GPOPS com BVP5C?','Sim', 'Não','Não');   
                       
    if strcmp(ButtonName6,'Sim')==1

    td = solution.time;
    ramseyBVP5C; 
    t = linspace(t0,tf,length(td)); 
    y = deval(sol,t);
    kbvp = y(1,:)';
    cbvp = y(2,:)';
    lambdabvp = costate(cbvp,rho,sigma,n,gx,t);
    
    figure
    subplot(2,2,1);
    plot(td, k,'b', t, kbvp,'r','LineWidth',1.5);
    set(gca,'FontName','Times','FontSize',16);
    xlabel('t','Interpreter','latex');
    ylabel('$k_a(t)$','Interpreter','latex');
    legend('k_a(t) GPOPS', 'k_a(t) BVP5C');
    title('GPOPS vs BVP5C - Capital');
    grid on;

    subplot(2,2,2);
    plot(td, c,'b' ,t, cbvp, 'r','LineWidth',1.5);
    set(gca,'FontName','Times','FontSize',16);
    xlabel('t','Interpreter','latex');
    ylabel('$c_a(t)$','Interpreter','latex');
    legend('c_a(t) GPOPS', 'c_a(t) BVP5C');
    title('GPOPS vs BVP5C - Consumo');
    grid on;
    assignin('base','kBVP5C',kbvp);
    assignin('base','cBVP5C',cbvp);
    
    subplot(2,2,3)
    plot(td,-solution.costate,'b',t, lambdabvp,'r','LineWidth',1.5);
    set(gca,'FontName','Times','FontSize',16);
    xlabel('$t$','Interpreter','latex');
    ylabel('$\lambda(t)$','Interpreter','latex');
    title('GPOPS vs BVP5C - Multiplicador');
    legend('Co-estado GPOPS','Co-estado BVP5C');
    grid on;
    assignin('base','lambdaBVP5C', lambdabvp);
        
        if strcmp(ButtonName5,'Sim') == 1 && alpha == 1/sigma
            subplot(2,2,1);
            hold on;
            plot(td, ksol,'k','LineWidth',1.5);
            legend('k_a(t) GPOPS', 'k_a(t) BVP5C','k_a(t) sol. analítica')
            hold off
            subplot(2,2,2);
            hold on
            plot(td, csol,'k','LineWidth',1.5);
            legend('c_a(t) GPOPS', 'c_a(t) BVP5C','c_a(t) sol. analítica')
            hold off
            subplot(2,2,3);
            hold on
            plot(td, lambdasol,'k','LineWidth',1.5);
            legend('\lambda GPOPS', '\lambda BVP5C','\lambda sol. analítica')
            hold off            
        else
        end

    else
    end
else
end
%---------------------------------------------------------------%

function ramseyGPOPS      
%Método Directo - GPOPS 
iphase = 1;

limits(iphase).time.min = [t0 tf];
limits(iphase).time.max = [t0 tf];
% Exclusão de trajectória explosiva que viola eq. de Euler k(T) = 0 dc/dt = +inf.
limits(iphase).state.min = [k0 0 kss];
limits(iphase).state.max = [k0 kss kss];
limits(iphase).control.min    =  0;
limits(iphase).control.max    =  100;
limits(iphase).parameter.min  = [];
limits(iphase).parameter.max  = [];
limits(iphase).path.min       = [];
limits(iphase).path.max       = [];
limits(iphase).event.min      = [];
limits(iphase).event.max      = [];
limits(iphase).duration.min   = [];
limits(iphase).duration.max   = [];

guess(iphase).time            = [t0; tf];
guess(iphase).state           = [k0; kss];
guess(iphase).control         = [c0; css];
guess(iphase).parameter       = [];

setup.name  = 'ramseymain-Problem';
setup.funcs.cost = 'ramseygpopsCost';
setup.funcs.dae = 'ramseygpopsDae';
setup.limits = limits;
setup.guess = guess;
setup.linkages = [];
setup.derivatives = 'finite-difference';
setup.autoscale = 'off';
setup.mesh.tolerance = 1e-3;
setup.mesh.iteration = 20;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 12;
setup.printoff = 1;

[output] = gpops(setup);
solution = output.solution;
solutionPlot = output.solutionPlot;
assignin('base','solutionGPOPS',solution);
assignin('base','solutionPlotGPOPS',solutionPlot);

end

function ramseyBVP5C
        
%Método Indirecto - bvp5c
NMax = 60;
RelTol = 1e-7;
AbsTol = 1e-4;
guessbvp = [k0 c0];     
options = bvpset('FJacobian',@J,'RelTol',RelTol,...
    'AbsTol',AbsTol','NMax',NMax);
solinit = bvpinit(linspace(t0,tf,10),guessbvp);
sol = bvp5c(@eqdinamicas,@bcfun,solinit,options);
end

function lambda = costate(c,rho,sigma,n,gx,t)
    %Coestado deduzida a partir do cálculo de variações (método indireto)
    l = length(c);
    lambda = zeros(l,1);
    for i = 1:l
        lambda1 = exp(gx.*t(i).*(1-1/sigma)-(rho-n).*t(i)).*c(i).^(-1/sigma);
        lambda(i) = lambda1;
    end
end

function ydot = eqdinamicas(~,y)
%equaçoes diferenciais com y(1)=k e y(2)=c
ydot = zeros(2,1);
ydot(1) = y(1).^alpha-(delta+n+gx).*y(1)-y(2); %dk/dt
ydot(2) = y(2).*sigma.*(alpha.*y(1).^(alpha-1)-delta-rho-gx/sigma); %dc/dt
end

function BGPka = kgr(func,k,c)  
%taxa de crescimento ( ou ka)   
l = length(k);
kdot = zeros(l,1);
BGPka = zeros(l,1);

    for m = 1:l

        ydot = feval(func,0,[k(m); c(m)]);
        kdot(m) = ydot(1);
        BGPka(m) = kdot(m)./k(m);
    end
end

function BGPca = cgr(func,k,c)
   %taxa de crescimento c    
z = length(c);
cdot = zeros(z,1);
BGPca = zeros(z,1);

    for i = 1:z
        ydot = feval(func,0,[k(i); c(i)]);
        cdot(i) = ydot(2);
        BGPca(i) = cdot(i)./c(i);
    end
end

function campovectorial(func,kval,cval)

n1 = length(kval);
n2 = length(cval);
%yp1 e yp2 são os vectores associados a yval e rval, respectivamente.
yp1 = zeros(n2,n1);   %y1val e y2val têm que ser igualmente espaçados
yp2 = zeros(n2,n1);
for i = 1:n1
  for j = 1:n2
     %yp1 e yp2 são os vectores associados a yval e rval, respectivamente. 
    ypv = feval(func,0,[kval(i);cval(j)]);        
    yp1(j,i) = ypv(1);
    yp2(j,i) = ypv(2);
  end
end
%Representação gráfica; Desenha setas nos vectores de velocidade
quiver(kval,cval,yp1,yp2,'r','LineWidth',1.4);       
axis tight;

end

function res = bcfun(ya,yb)
%para o BVP5C
  res = [ ya(1) - k0
          yb(1) - kss];  
end

function j = J(~,y)

%Matriz jacobiana do sistema em @eqdinamicas
j = zeros(2,2);
j(1,1) = alpha.*y(1,:).^(alpha-1)-(delta+n+gx);
j(1,2) = -1;
j(2,1) = y(2,:).*sigma.*(alpha.*(alpha-1)).*y(1,:).^(alpha-2);
j(2,2) = sigma.*(alpha.*y(1,:).^(alpha-1)-delta-rho-gx/sigma);
end

function yka = gdp(k)
%yka = produto por unidade eficiente de trabalho
l = length(k);
yka = zeros(l,1);
    for i = 1:l
    gdp = k(i).^alpha;
    yka(i) = gdp;
    end

end

function s = savings(yka, c)
%poupança por unidade eficiente de trabalho
l = length(c);
s = zeros(l,1);
    for i = 1:l
        savings = yka(i) - c(i);
    s(i) = savings;
    end
end

end