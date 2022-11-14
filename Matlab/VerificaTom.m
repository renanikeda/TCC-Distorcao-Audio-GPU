clear all;
close all;

indice2=2;
for A=[0.1 0.3 0.5 0.7 0.9 1]
indice=1;
for grau=[1/10 1/5 1/3 1/2 1 2 5 10]

fs=100;
fa=2*(2*fs)*3*3*3;
N=(3)*round(fa/fs)
omegao=fs/fa;
n=0:N-1;
x=A*cos(2*pi*omegao*n)+A*cos(4*pi*omegao*n);

y=atan(grau*x)/atan(grau);
z=atan(grau*y)/atan(grau);
h=atan(grau*z)/atan(grau);
%  
eA=1.2;
nt=1000*n/fa; % ms
vetortempo=[0 max(nt) -eA eA];
figure(1)
subplot(241); plot(nt,x); title('x(n)')
axis(vetortempo); grid; xlabel('ms')
subplot(242); plot(nt,y); title('y(n)')
axis(vetortempo); grid; xlabel('ms')
subplot(243); plot(nt,z); title('z(n)')
axis(vetortempo); grid; xlabel('ms')
subplot(244); plot(nt,h); title('h(n)')
axis(vetortempo); grid; xlabel('ms')
% 
passo=1/1024;
omega=(0:passo:(1-passo))*2*pi;
K=length(omega);
K2=round(K/2);

% 
%w = window(@hamming,N)';
w = window(@hanning,N)';
%w = window(@blackman,N)';
%w=ones(1,N);

X=20*log10(abs(fft(x.*w,K)));
Y=20*log10(abs(fft(y.*w,K)));
Z=20*log10(abs(fft(z.*w,K)));
H=20*log10(abs(fft(h.*w,K)));
W=20*log10(abs(fft(w,K)));
OM=(fa/2)*[0:K2-1]/K2;
vetor=[0 3000 -100 50];
subplot(245); plot(OM,X(1:K2));hold;plot(OM,W(1:K2))
axis(vetor); 
legend('X','W');
title('TFTD {x(n)}');grid
subplot(246); plot(OM,Y(1:K2));hold;plot(OM,W(1:K2))
axis(vetor); 
legend('Y','W');
title('TFTD {y(n)}');grid
subplot(247); plot(OM,Z(1:K2));hold;plot(OM,W(1:K2))
axis(vetor); 
legend('Z','W');
title('TFTD {z(n)}');grid
subplot(248); plot(OM,H(1:K2));hold;plot(OM,W(1:K2))
axis(vetor); 
legend('H','W');
title('TFTD {h(n)}');grid

figure(indice2)
subplot(4,2,indice)
[0 K2-1]
plot(OM,X(1:K2), 'LineWidth',2);
hold on
plot(OM,Y(1:K2), 'LineWidth',2);
plot(OM,Z(1:K2), 'LineWidth',2);
plot(OM,H(1:K2), 'LineWidth',2);
hold off
legend('1 estágio','2 estágio','3 estágio','4 estágio')
title(['Espectro de atan(Acos(alfa*x)) para ','A = ',num2str(A),' grau = ',num2str(grau)]);
grid
xlabel('Hz') 
ylabel('Módulo[dB]')
axis([0 3000 -20 50])

indice=indice+1;
end
indice2=indice2+1;
end