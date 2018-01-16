close all
clear;

load('fcno01fz.mat');
%% Q-1 :

fe=8000;  
N=1025;
sigma=0.25;
bruit=sigma*randn(1,N);

%représentation théorique de la fonction d'autocorrelation.

x=(-floor(N/2):1:floor(N/2));
y=zeros(1,N);
y(floor(N/2)+1)=sigma^2;
% figure
% plot(x,y)

%représentation estimée

% biased_br=xcorr(bruit,'biased');
% unbiased_br=xcorr(bruit,'unbiased');
% N1=length(biased_br);
% x1=(-floor(N1/2):1:floor(N1/2));
% figure
% subplot(2,1,1)
% plot(x1,biased_br,'r')
% subplot(2,1,2)
% plot(x1,unbiased_br,'g')


%spectre de puissance d'une réalisation & densité spectrale de puissance du
%bruit.

fbruit=fft(bruit);
DSPy=fbruit.*conj(fbruit);                                                 
SP=((abs(fbruit)).^2)./N;
f=(-1/2:1/N:1/2-1/N);
 
% figure
% subplot(2,1,1)
% plot(f,DSPy,'r')
% subplot(2,1,2)
% plot(f,SP,'b')

%% préliminaire2

s = fcno01fz';
RSB = 10;
len = length(s);
s_noise  = noise_sig( s, RSB );   % On ajoute au signal un BBC de RSB=5.

t = 0:1/fe:len/fe-1/fe;
figure
subplot(2,1,1)
plot(t, s);
xlim([0 6.5 ]);
xlabel('temps (s)');
ylabel('s');
title('representation temporelle du signal original');
subplot(2,1,2)

spectrogram(s,256,0,256,fe,'yaxis')
axis tight;
xlim([0 6.5 ]);



figure
subplot(2,1,1)
plot(t,s_noise)
xlim([0 6.5 ]);
xlabel('temps (s)');
ylabel('s bruité');
title('representation temporelle du signal bruite');

subplot(2,1,2)

spectrogram(s_noise,256,0,256,fe,'yaxis')
xlim([0 6.5 ]);

%% Procedure d'addition-recouvrement :

s_noise1=noise_sig(s, 5);

%découpage
taux_rec = 0.5;  % Taux de recouvrement est de 50%.
len_trame=256;
nbr_trames=floor((len-len_trame)/(len_trame-taux_rec*len_trame))+1;
trame_matrix = zeros(nbr_trames, len_trame);
win_han=hanning(256)';
trame_matrix(1,:) = (1:len_trame).*win_han;

f_1 = len_trame-taux_rec*len_trame;  % Chaque fois on decale d'un moitié de trame.
f_2 = len_trame-taux_rec*len_trame;

for i=2:nbr_trames
    
        trame_matrix(i,:) = s(f_2+1:f_2+len_trame);
        trame_matrix(i,:)=trame_matrix(i,:).*win_han;
        f_2 = f_2+f_1;
    
end


%% Recupération

%moitié de la premiere trame et la derniere trame
trame_matrix_recup=zeros(nbr_trames, len_trame);
trame_matrix_recup(1,(1:f_1))=trame_matrix(1,(1:f_1));
trame_matrix_recup(nbr_trames,(f_1+1:end))=trame_matrix(nbr_trames,(f_1+1:end));

%pour le reste moyennage
s_recup=zeros(1,len);
s_recup(1:f_1)=trame_matrix_recup(1,(1:f_1));
s_recup(end-f_1+1:end)=trame_matrix(nbr_trames,(f_1+1:end));
f_2 = f_1;

for i=1:(len-len_trame)/f_1
    
   s_recup(f_2+1:(len_trame-f_1)+f_2) = (trame_matrix(i,(f_1+1:end))+trame_matrix(i+1,(1:f_1)))./(win_han(1:f_1)+win_han(f_1+1:end)); 
   f_2 = f_2 + f_1;   
   
end

% Comparaison entre signal recuperer et signal original
len = length(s);
t = 0:1/fe:len/fe-1/fe;
figure
plot(t,s_recup);
xlim([ 0 6.5 ]);
xlabel('temps (s)');
ylabel('Signal reconstruit');

figure
plot(t,s)
xlim([ 0 6.5 ]);
xlabel('temps (s)');
ylabel('Signal originale');

e = s-s_recup;        % Le calcule du taux d'erreur.
%b=sum(abs(e(f_1+1:end-f_1)));
%error=b/length(e)













