%Interpolação por filtro IIR sem blocos

clear all;
close all;

%Definições Iniciais
L=30;
N_sinal=48000*5;
N_fft=N_sinal*L;
N_pulso=N_sinal*L;
grau_distorcao=0.2;
w = window(@hanning,N_sinal)';
w_interpolado = window(@hanning,N_sinal*L)';
N_estagios=5;

%Filtro de Interpolacao
Rp_interpol=1;
Rs_interpol=90;
Wp_interpol=2*15000/(L*48000);
Ws_interpol=1/L;
[N_interpol,Wp_interpol] = ellipord(Wp_interpol,Ws_interpol,Rp_interpol,Rs_interpol);
[b_interpol,a_interpol] = ellip(N_interpol, Rp_interpol, Rs_interpol, Wp_interpol);
[H_interpol,W]=freqz(b_interpol,a_interpol,N_fft/2+1);

%Filtro de Decimacao
N_decimacao=N_interpol;
a_decimacao=a_interpol;
b_decimacao=b_interpol;

%Filtro anti Aliasing
Rp_aliasing=0.05/1000;
Rs_aliasing=50;
Wp_aliasing=2*15000/(L*48000);
Ws_aliasing=1/L;
[N_aliasing,Wp_aliasing] = ellipord(Wp_aliasing,Ws_aliasing,Rp_aliasing,Rs_aliasing);
[b_aliasing,a_aliasing] = ellip(N_aliasing, Rp_aliasing, Rs_aliasing, Wp_aliasing);
[H_aliasing,W_aliasing]=freqz(b_aliasing,a_aliasing,N_fft/2+1);

%Leitura do sinal e espectro
x=audioread('sweep_48khz_5s.wav');
x=x(1:N_sinal); %Audio Gravado
x_espectro=fft(x'.*w,N_fft);

%Começo
entrada_interpol=zeros(1,N_interpol+1);
saida_interpol=zeros(1,N_interpol);
entrada_aliasing=zeros(N_estagios,N_aliasing+1);
saida_aliasing=zeros(N_estagios,N_aliasing);
entrada_decimacao=zeros(1,N_decimacao+1);
saida_decimacao=zeros(1,N_decimacao);

x_interpolado_filtrado_total=zeros(1,N_sinal*L);
x_distorcido_total=zeros(1,N_sinal*L);
x_superamostrado=zeros(1,L*N_sinal);
x_distorcido_filtrado_total=zeros(1,L*N_sinal);
x_pre_decimacao_total=zeros(1,L*N_sinal);
x_decimado_total=zeros(1,N_sinal);
x_superamostrado=repelem(x,L);
a_interpol=a_interpol(2:N_interpol+1);
a_aliasing=a_aliasing(2:N_aliasing+1);
a_decimacao=a_decimacao(2:N_decimacao+1);
aux=1;

%Salvar os coeficientes dos filtros
arqcoef(a_interpol,"a_interpol_30.txt");
arqcoef(b_interpol,"b_interpol_30.txt");
arqcoef(a_decimacao,"a_decimacao_30.txt");
arqcoef(b_decimacao,"b_decimacao_30.txt");
arqcoef(a_aliasing,"a_aliasing_30.txt");
arqcoef(b_aliasing,"b_aliasing_30.txt");

for i=1:L*N_sinal

%Interpolação
for j=N_interpol+1:-1:2
    entrada_interpol(j)=entrada_interpol(j-1);
end
entrada_interpol(1)=x_superamostrado(i);

y_interpol=sum(b_interpol.*entrada_interpol)-sum(a_interpol.*saida_interpol);

for j=N_interpol:-1:2
    saida_interpol(j)=saida_interpol(j-1);
end
saida_interpol(1)=y_interpol;

if y_interpol>1
     y_interpol=1;
end
if y_interpol<-1
     y_interpol=-1;   
end

x_interpolado_filtrado_total(i)=y_interpol;


%Distorcao com filtragem
for distorcoes=1:N_estagios
    
if distorcoes==1
    x_distorcido=atan(y_interpol.*grau_distorcao)/atan(grau_distorcao);
else
    x_distorcido=atan(y_aliasing.*grau_distorcao)/atan(grau_distorcao);
end

for j=N_aliasing+1:-1:2
    entrada_aliasing(distorcoes,j)=entrada_aliasing(distorcoes,j-1);
end
entrada_aliasing(distorcoes,1)=x_distorcido;

y_aliasing=sum(b_aliasing.*entrada_aliasing(distorcoes,:))-sum(a_aliasing.*saida_aliasing(distorcoes,:));

for j=N_aliasing:-1:2
    saida_aliasing(distorcoes,j)=saida_aliasing(distorcoes,j-1);
end
saida_aliasing(distorcoes,1)=y_aliasing;

end

x_distorcido_filtrado_total(i)=y_aliasing;

%Decimacao
for j=N_decimacao+1:-1:2
    entrada_decimacao(j)=entrada_decimacao(j-1);
end

%entrada_decimacao(1)=y_interpol;
entrada_decimacao(1)=y_aliasing;

y_decimacao=sum(b_decimacao.*entrada_decimacao)-sum(a_decimacao.*saida_decimacao);


for j=N_decimacao:-1:2
    saida_decimacao(j)=saida_decimacao(j-1);
end
saida_decimacao(1)=y_decimacao;
x_pre_decimacao_total(i)=y_decimacao;

%Decimacao
if mod(i-1,L)==0
    x_decimado_total(aux)=y_decimacao;
    if y_decimacao>1
        x_decimado_total(aux)=1;
    else
        if y_decimacao<-1
            x_decimado_total(aux)=-1;
        end
    end
    aux=aux+1;
end

end

% %Plots
figure(1)
subplot(3,1,1)
plot(1:N_sinal,x);
subplot(3,1,2)
plot(1:N_sinal*L,x_interpolado_filtrado_total);
subplot(3,1,3)
plot(1:N_sinal,x_decimado_total);

x_interpolado_filtrado_total_espectro=fft(x_interpolado_filtrado_total.*w_interpolado,N_fft);
x_decimado_total_espectro=fft(x_decimado_total.*w,N_fft);

figure(2)
subplot(2,1,1)
plot((0:(N_fft/2))/(N_fft/2),20*log10(abs(x_espectro(1:(N_fft/2)+1))));
hold
plot((0:(N_fft/2))/(N_fft/2),20*log10(abs(x_decimado_total_espectro(1:(N_fft/2)+1))));
legend('Original','Decimado Total')
subplot(2,1,2)
plot((0:(N_fft/2))/(N_fft/2),20*log10(abs(x_interpolado_filtrado_total_espectro(1:(N_fft/2)+1))));


%Audição
% sound(x,44100)
% pause
%sound(x_decimado_total,44100)