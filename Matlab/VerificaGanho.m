close all;
clear all;

%Verifica Ganho 

N_pontos=100000;
n=(-N_pontos:1:N_pontos)/N_pontos;
N_distorcoes=2000;
N_distorcoes=N_distorcoes-1;
limite_erro=10^-1;

%Retangular
retangular=zeros(1,length(n));
for i=1:length(n);
    if n(i)>0
        retangular(i)=1;
    else
        retangular(i)=-1;
    end
end

figure(1)
grau_distorcao=0.05;
distorcao=atan(n.*grau_distorcao)/atan(grau_distorcao);
plot(n,retangular);
hold
plot(n,distorcao);


for grau_distorcao=1.0:0.5:5.0
    distorcao=atan(n.*grau_distorcao)/atan(grau_distorcao);
    plot(n,distorcao);
end

grid;
title('Distorcão atan(grau*x)/atan(grau)')
xlabel('x') 
ylabel('atan(grau*x)/atan(grau)')
legend('sign(x)','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0','4.5','5.0');

%Equivalente Composto
for i=1:N_distorcoes
    distorcao=atan(distorcao.*grau_distorcao)/atan(grau_distorcao);
end

figure(2)
plot(n,retangular);
hold
grau_distorcao=0.05;
distorcao=atan(n.*grau_distorcao)/atan(grau_distorcao);
for i=1:N_distorcoes
    distorcao=atan(distorcao.*grau_distorcao)/atan(grau_distorcao);
end
plot(n,distorcao);




for grau_distorcao=0.2:0.1:0.5
    distorcao=atan(n.*grau_distorcao)/atan(grau_distorcao);
    for i=1:N_distorcoes
        distorcao=atan(distorcao.*grau_distorcao)/atan(grau_distorcao);
    end
    plot(n,distorcao);
end

grid;
title(['Distorcão atan(grau*x)/atan(grau) composto para ', num2str(N_distorcoes+1),' estágios']);
xlabel('x') 
ylabel(['atan(grau*x)/atan(grau)composto para ', num2str(N_distorcoes+1),' estágios']);
legend('sign(x)','0.1','0.2','0.3','0.4','0.5');


for j=1:20
    grau_distorcao=1/j;
    distorcao=n;
    erro(j)=1;
    N_necessario(j)=0;
    while erro(j)>limite_erro
        distorcao=atan(distorcao.*grau_distorcao)/atan(grau_distorcao);
        erro(j)=sum((distorcao-retangular).^2)/length(n);
        N_necessario(j)=N_necessario(j)+1;
    end
end

figure(3)
semilogy(1./(1:20),N_necessario);
grid;
title(['Numero de estágios necessarios para EQM <= ',num2str(limite_erro)]);
xlabel('grau de distorção') 
ylabel('Numero de estágios necessarios');

