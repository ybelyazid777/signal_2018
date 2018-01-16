%% Procedure d'addition-recouvrement :

close all

load('fcno01fz.mat');
s_noise1=noise_sig(s, 5);
s = fcno01fz';
f_e=8000;

%découpage
taux_rec = 0.5;
len_trame=256;
nbr_trames=floor((len-len_trame)/(len_trame-taux_rec*len_trame))+1;
trame_matrix = zeros(nbr_trames, len_trame);
trame_matrix_b = zeros(nbr_trames, len_trame);
win_han=hanning(256)';
trame_matrix(1,:) = (1:len_trame).*win_han;
trame_matrix_b(1,:) = (1:len_trame).*win_han;
f_1 = len_trame-taux_rec*len_trame;
f_2 = len_trame-taux_rec*len_trame;

for i=2:nbr_trames
    
        trame_matrix(i,:) = s(f_2+1:f_2+len_trame); 
        trame_matrix(i,:)=trame_matrix(i,:).*win_han;
        trame_matrix_b(i,:) = s_noise1(f_2+1:f_2+len_trame);
        trame_matrix_b(i,:)=trame_matrix_b(i,:).*win_han;
        f_2 = f_2+f_1;
    
end


%% estimation de la variance du bruit

% figure, plot(s_noise1)
h=s_noise1(end-256+1:end);
DSP_b=0;

for i=1:40
    
    h = s_noise1(end-i*256+1:end-(i-1)*256);
    DSP_b = DSP_b+var(h);
    
end

DSP_b = DSP_b/40;

%reauhaussement

fft_trame_b=zeros(nbr_trames, len_trame);
phase_trame_b=zeros(nbr_trames, len_trame);
sp_trame_b=zeros(nbr_trames, len_trame);
sp_trame_rehaus=zeros(nbr_trames, len_trame);
module_trame_rehaus=zeros(nbr_trames, len_trame);
fft_trame_rehaus=zeros(nbr_trames, len_trame);
s_rehaus=zeros(nbr_trames, len_trame);


for i=1:nbr_trames
    fft_trame_b(i,:)=fft(trame_matrix_b(i,:));
    phase_trame_b(i,:)=angle(fft_trame_b(i,:));
    sp_trame_b(i,:)=(abs(fft_trame_b(i,:).^2)/length(fft_trame_b(i,:)));
    sp_trame_rehaus(i,:)=sp_trame_b(i,:)-DSP_b;
    sp_trame_rehaus(sp_trame_rehaus<0)=0;
    module_trame_rehaus(i,:)=sqrt(sp_trame_rehaus(i,:).*length(sp_trame_rehaus(i,:)));
    fft_trame_rehaus(i,:)=module_trame_rehaus(i,:).*(exp(j.*phase_trame_b(i,:)));
    s_rehaus(i,:)=ifft(fft_trame_rehaus(i,:));
    
%     spbruite=(abs(S_noise1).^2)/length(S_noise1);
%     sprehaus=spbruite-DSP_b;
%     sprehaus(sprehaus<0)=0;
% 
%     module_rehaus=sqrt(sprehaus*length(sprehaus));
% 
%     fft_rehaus=module_rehaus.*(exp(j.*v_phase));
% 
%     signal_rehaus=ifft(fft_rehaus,'symmetric');

end
    
%recupération

%moitié de la premiere trame et la derniere trame
trame_matrix_recup=zeros(nbr_trames, len_trame);
trame_matrix_recup(1,(1:f_1))=s_rehaus(1,(1:f_1));
trame_matrix_recup(nbr_trames,(f_1+1:end))=s_rehaus(nbr_trames,(f_1+1:end));

%pour le reste moyennage
s_recup=zeros(1,len);
s_recup(1:f_1)=trame_matrix_recup(1,(1:f_1));
s_recup(end-f_1+1:end)=s_rehaus(nbr_trames,(f_1+1:end));
f_2 = f_1;

for i=1:(len-len_trame)/f_1
    
   s_recup(f_2+1:(len_trame-f_1)+f_2) = (s_rehaus(i,(f_1+1:end))+s_rehaus(i+1,(1:f_1)))./(win_han(1:f_1)+win_han(f_1+1:end)); 
   f_2 = f_2 + f_1;    
end

len = length(s);
t = 0:1/fe:len/fe-1/fe;
figure
plot(t,real(s_recup));
xlim([ 0 6.5 ]);
xlabel('temps (s)');
ylabel('Signal recupere');


figure
plot(t,s)
xlim([ 0 6.5 ]);
xlabel('temps (s)');
ylabel('Signal original');


figure
plot(t,s_noise1);
xlim([ 0 6.5 ]);
xlabel('temps (s)');
ylabel('Signal bruite');

t_trame = 0:1/f_e:2*f_1/f_e-1/f_e

figure

plot(t_trame,trame_matrix(20,:),'g');
hold on
plot(t_trame,trame_matrix_b(20,:),'r');
plot(t_trame,s_rehaus(20,:),'b');
hold off

legend('trame originale','trame bruitée','trame rehaussée(x)' );
title('representation temporelle : trame originale, bruitée et rehaussée');
xlabel('temps (s)');
xlim([ 0 0.03 ]);




