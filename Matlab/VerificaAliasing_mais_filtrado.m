%Interpolação por filtro IIR sem blocos

clear all;
close all;

%Definições Iniciais
aux=1;
for L=2:2:30
posicao=1;
for grau_distorcao=[1/10 1/5 1/3 1/2 1]

N_sinal=48000*0.1;
N_fft=N_sinal*L;
w = window(@hanning,N_sinal)';
w_interpolado = window(@hanning,N_sinal*L)';


%Filtro
Rp=0.01;
Rs=33;
Wp=2*20000/(L*48000);
Ws=1/L;
[N,Wp] = ellipord(Wp,Ws,Rp,Rs);
[b,a] = ellip(N, Rp, Rs, Wp);
% b_interpol=L.*b;
b_interpol=b;
[H,W]=freqz(b,a,N_fft/2+1);

%Leitura do sinal e espectro
x=audioread('sweep_48khz_100ms.wav');
x=x(1:N_sinal); %Audio Gravado
x_espectro=fft(x'.*w,N_fft);

%Começo
x_interpolado_filtrado_total=zeros(1,L*N_sinal);
entrada_primeiro=zeros(1,N+1);
saida_primeiro=zeros(1,N);
entrada_segundo=zeros(1,N+1);
saida_segundo=zeros(1,N);
entrada_terceiro=zeros(1,N+1);
saida_terceiro=zeros(1,N);
x_superamostrado=zeros(1,L*N_sinal);
% x_superamostrado(1:L:end)=x;
x_superamostrado=repelem(x,L);
x_interpolado_filtrado_distorcido_total=[];
a=a(2:N+1);

for i=1:L*N_sinal

%Primeira Filtragem
for j=N+1:-1:2
    entrada_primeiro(j)=entrada_primeiro(j-1);
end
entrada_primeiro(1)=x_superamostrado(i);
    
y_primeiro=sum(b_interpol.*entrada_primeiro)-sum(a.*saida_primeiro);

for j=N:-1:2
    saida_primeiro(j)=saida_primeiro(j-1);
end
saida_primeiro(1)=y_primeiro;

%Segunda Filtragem
for j=N+1:-1:2
    entrada_segundo(j)=entrada_segundo(j-1);
end
entrada_segundo(1)=y_primeiro;
    
y_segundo=sum(b_interpol.*entrada_segundo)-sum(a.*saida_segundo);

for j=N:-1:2
    saida_segundo(j)=saida_segundo(j-1);
end
saida_segundo(1)=y_segundo;

%Terceira Filtragem
for j=N+1:-1:2
    entrada_terceiro(j)=entrada_terceiro(j-1);
end
entrada_terceiro(1)=y_segundo;
    
y_terceiro=sum(b_interpol.*entrada_terceiro)-sum(a.*saida_terceiro);

for j=N:-1:2
    saida_terceiro(j)=saida_terceiro(j-1);
end
saida_terceiro(1)=y_terceiro;

if(y_terceiro>1)
    y_terceiro=1;
end

if(y_terceiro<-1)
    y_terceiro=-1;
end
    
x_interpolado_filtrado_total(i)=y_terceiro;

%Distorcao
x_interpolado_filtrado=atan(y_terceiro.*grau_distorcao)/atan(grau_distorcao);
x_interpolado_filtrado_distorcido_total=[x_interpolado_filtrado_distorcido_total x_interpolado_filtrado];

end

x_interpolado_filtrado_total_espectro=fft(x_interpolado_filtrado_total.*w_interpolado,N_fft);
x_interpolado_filtrado_distorcido_total_espectro=fft(x_interpolado_filtrado_distorcido_total.*w_interpolado,N_fft);
x_aliasing=x_interpolado_filtrado_distorcido_total_espectro-x_interpolado_filtrado_total_espectro;


maximo_passagem(posicao,aux)=max(20*log10(abs(x_interpolado_filtrado_total_espectro(1:(N_fft/(2*L))+1))));
maximo_passagem_com_distorcao(posicao,aux)=max(20*log10(abs(x_interpolado_filtrado_distorcido_total_espectro(1:(N_fft/(2*L))+1))));
maximo_pi(posicao,aux)=max(20*log10(abs(x_aliasing(round(0.7*(N_fft/(2)+1)):(N_fft/(2)+1)))));
aliasing(posicao,aux)=maximo_passagem_com_distorcao(posicao,aux)-maximo_pi(posicao,aux);

figure(aux)
if posicao==1
plot((0:(N_fft/2))/(N_fft/2),20*log10(abs(x_interpolado_filtrado_total_espectro(1:(N_fft/2)+1)))); 
hold
plot((0:(N_fft/2))/(N_fft/2),20*log10(abs(H(1:(N_fft/2)+1))));
plot((0:(N_fft/2))/(N_fft/2),20*log10(abs(x_interpolado_filtrado_distorcido_total_espectro(1:(N_fft/2)+1))));
else
plot((0:(N_fft/2))/(N_fft/2),20*log10(abs(x_interpolado_filtrado_distorcido_total_espectro(1:(N_fft/2)+1)))); 
end

if posicao==5
grid;
title('Espectro interpolado e distorcido do sweep')
xlabel('frequência discreta normalizada') 
ylabel('Módulo [dB]')
axis([0 1 -200 100]);
legend('Interpolado','Filtro','Grau distocao = 0.1','Grau distocao = 0.2','Grau distocao = 0.3','Grau distocao = 0.5','Grau distocao = 1.0')
end

posicao=posicao+1;
end

aux=aux+1;
L
end



%Plots
figure(aux)
plot(2:2:30,maximo_pi(1,1:15));
hold
for indice=2:5
    plot(2:2:30,maximo_pi(indice,1:15));
end
legend('Grau distocao = 0.1','Grau distocao = 0.2','Grau distocao = 0.33','Grau distocao = 0.5','Grau distocao = 1.0')
title('Espectro em pi para diferentes L e graus de distorcao')
xlabel('L') 
ylabel('Modulo do espectro')

figure(aux+1)
plot(2:2:30,maximo_passagem_com_distorcao(1,1:15));
hold
for indice=2:5
    plot(2:2:30,maximo_passagem_com_distorcao(indice,1:15));
end
legend('Grau distocao = 0.1','Grau distocao = 0.2','Grau distocao = 0.33','Grau distocao = 0.5','Grau distocao = 1.0')
title('Espectro na faixa de passagem distorcido para diferentes L e graus de distorcao')
xlabel('L') 
ylabel('Modulo do espectro')

figure(aux+2)
plot(2:2:30,aliasing(1,1:15));
hold
for indice=2:5
    plot(2:2:30,aliasing(indice,1:15));
end
legend('Grau distocao = 0.1','Grau distocao = 0.2','Grau distocao = 0.33','Grau distocao = 0.5','Grau distocao = 1.0')
title('Diferença entre espectro em pi e máximo da faixa de passagem para vários L e graus de distorcao')
xlabel('L') 
ylabel('Diferença em dB')


