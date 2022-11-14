clear all;
close all;

for A=1:-0.1:0.1
%Verifica formulas de cosseno
a=1;
% A=0.5;
n=(-1000:1000)/2;

%Matriz de Harmonicos
N = 1000;
M = zeros(N,N+1);

M(1,N) = 1;
M(2,N-1) = 1;

for i = 3:N
    for j = N:-1:2
        M(i,j) = M(i-1,j-1) + M(i-1, j+1);
        
        if (j == N-1)
            if (mod(i,2) == 0)
                M(i,j) =  M(i-1,j-1) + M(i-1, j+1) * 2;
            end
        end
    end
end

M(N,1) = 1;
M = M(:, 1:N);
for i = 2:N
    M(i,:) = A^(i-1).*M(i,:)./(2^(i-2));
end

%Gerado
real=(A.*cos(2*n*pi/11)).^(N-1);
formula=zeros(1,2001);
for k=1:N
    formula=formula+(M(N,k)*(cos((N-k)*2*n*pi/11)));
end

figure(1)
plot(n,real)
hold
plot(n,formula)
legend('Cos^N-1','Formula')

%Coeficientes de Taylor
Taylor = M;

indice = 1;

for i = 2:2:N
    Taylor(i,:) = a^(i-1).*Taylor(i,:)./(indice);
    
    if (indice > 0)
        indice = (indice+2)*(-1);
    else
        indice = (indice-2)*(-1);
    end
end

Coeficientes_harmonicas(round(A*10),:)=zeros(1,N-1);
indice=1;
for k=N-1:-2:1
    Coeficientes_harmonicas(round(A*10),indice)=sum(Taylor(:,k));
    indice=indice+1;
end

Coeficientes_harmonicas_dB(round(A*10),:)=20*log10(abs(Coeficientes_harmonicas(round(A*10),:)));

taylor_harm=zeros(1,2001);
for i=1:length(Coeficientes_harmonicas(round(A*10),:))
    taylor_harm=taylor_harm+Coeficientes_harmonicas(round(A*10),i).*(cos((2*i-1)*2*n*pi/9));
end

real=atan(a*A*cos(2*n*pi/9));
taylor_pot=zeros(1,2001);
x=A*cos(2*n*pi/9);
for order=0:499
taylor_pot=taylor_pot+(((-1)^(2+order))*(x.^(2*order+1)))/((2*order+1));
end

figure(2)
plot(n,taylor_harm)
hold
plot(n,real)
plot(n,taylor_pot)
legend('Taylor Harmonicos','Atan','Taylor Pot')


end

figure(3)
semilogx(2.*(0:length(Coeficientes_harmonicas_dB(10,:))-1)+1,Coeficientes_harmonicas_dB(9,:))
hold
for i=8:-1:1
    semilogx(2.*(0:length(Coeficientes_harmonicas_dB(i,:))-1)+1,Coeficientes_harmonicas_dB(i,:))
end
grid;
legend('A=0.9','A=0.8','A=0.7','A=0.6','A=0.5','A=0.4','A=0.3','A=0.2','A=0.1')
title('Coeficientes Harmonicas')
xlabel('enésima harmonica') 
ylabel('Amplitude em dB')
axis([1 100 -200 0]);

% %Frequencias rebatidas
% frequencia_harmonica=zeros(1,(N-2)/2+1);
% frequencia_harmonica(1)=1/9;
% 
% for i=1:499
%     harmonico=(1/9)*(2*i+1);
%     if (harmonico>1 && rem(harmonico,2)>=1)
%         frequencia_harmonica(i+1)=2-rem(harmonico,2);
%     else
%         if (harmonico>1 && rem(harmonico,2)<1)
%             frequencia_harmonica(i+1)=rem(harmonico,2);
%         else
%             frequencia_harmonica(i+1)=harmonico;
%         end
%     end
% end
% 
% %Juntando frequencias rebatidas
% Coeficientes_harmonicas_rebatidos=zeros(1,5);
% 
% for i=1:5
%     for j=1:500
%         if abs(frequencia_harmonica(i)-frequencia_harmonica(j))<0.01
%             Coeficientes_harmonicas_rebatidos(i)=Coeficientes_harmonicas_rebatidos(i)+Coeficientes_harmonicas(j);
%         end
%     end
% end
% 
% sinal_rebatido=zeros(1,2001);
% for i=1:5
%     sinal_rebatido=sinal_rebatido+Coeficientes_harmonicas_rebatidos(i)*cos((2*i-1)*2*pi*n/9);
% end
% 
% N_fft=10000;
% sinal_rebatido_fft=fft(sinal_rebatido,N_fft);
% real_fft=fft(real,N_fft);
% 
% figure(4)
% plot((0:(N_fft/2))/(N_fft/2),20*log10(abs(real_fft(1:(N_fft/2)+1))));
% hold
% plot((0:(N_fft/2))/(N_fft/2),20*log10(abs(sinal_rebatido_fft(1:(N_fft/2)+1))));
% legend('Real','Rebatido simulado')