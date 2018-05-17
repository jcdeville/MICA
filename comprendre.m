clear all;
close all;

load('ecg_normal_1.mat')

%% parameters

%low pass filter
low_pass_num=[1 0 0 0 0 0 -2 0 0 0 0 0 1];
low_pass_den=[1 -2 1];

%high pass filter
high_pass_num=zeros(1,33);
high_pass_num(1,1)=-1;
high_pass_num(1,17)=32;
high_pass_num(1,18)=-32;
high_pass_num(1,33)=1;
high_pass_den=[1 -1];


after_lpf=filter(low_pass_num,low_pass_den,ecg);
after_hpf=filter(high_pass_num,high_pass_den,after_lpf);


% figure(1)
% hold on
% plot(ecg)
% plot(after_lpf)
% plot(after_hpf)
% legend('ecg','after_lpf','after_hpf');
% hold off

%filter derivative
deri=[1 2 0 -2 -1].*Fs/8;
deriv=filter(deri,1,after_hpf);

% squared signal
ssq=abs(deriv).^2;


%moving windows integration
N=0.10*Fs+1; % En moyenne la durée d'un SQR complex est 0.10 secondes

[nb_ligne,nb_col]=size(ssq);
for n=N:nb_col
    sum_1=0;
    for i=0:N-1
    sum_1=sum_1+ssq(1,n-i);
    end
    swni(1,n)=1/N*sum_1;
end

%% tresholding

% on établi un seuil pour lesquels le signal obtenu est égale à 0

thresholing = mean(swni);
for j=1:nb_col
    if (swni(1,j)<thresholing)
        swni(1,j)=0;
    end
end


%ajout de 1 � N pour avoir un retard entier car le retard est de (N-1)/2
%delai introduit par la window = 10, car en fait c'est une convolution avec une porte (la
%somme)


%% calcul du retard du filtre passe bande

%fvtool(low_pass_num,low_pass_den) %nous donne un retard de 5
%fvtool(high_pass_num,high_pass_den) %donne environ 16


%% retard induit par l'algorithme de pan an tompkins

retard_2 = zeros(1,50); % 50 est arbitraire tant qu'il soit supérieur au retard
retard_2(1,(N-1)/2+2+1+21) = 1; % retard windows + derivative + ban-pass filter
signal_de_base_retarde = conv(ecg,retard_2);

%% Recherche des R 

indice = 1;
indice_3 = 1;

while indice<nb_col
    max_local = 0;
    indice_2 = 0;
    indice_max = 0;
    if(swni(1,indice)~=0)
        while(swni(1,indice+indice_2)~=0)
            if(max_local<signal_de_base_retarde(1,indice+indice_2))
               indice_max = indice+indice_2;
               max_local = signal_de_base_retarde(1,indice_max);
            end
            indice_2 = indice_2 + 1;
        end
        if(indice_max~=0)
         max(1,indice_3) = indice_max;
         max(2,indice_3) = signal_de_base_retarde(1,indice_max);
         indice_3 = indice_3 + 1;

        end
        indice = indice+ indice_2;
    else 
        indice=indice+1;
    end
end

[nb_ligne_1,nb_col_2]=size(max);


%% Idée de damien : regardé l'amplitude moyenne d'un point R 
%% Idée de jc : regardé la distance moyenne entre deux points R


mean_max = mean(max(2,:));
var_max = var(max(2,:));
k = 1;
for i=1:nb_col_2
    if(max(2,i)>=mean_max-5*sqrt(var_max) & max(2,i)<=mean_max+5*sqrt(var_max))
        corrected_max(1,k)=max(1,i);
        corrected_max(2,k)=max(2,i);
        k=k+1;
    end
end
max = corrected_max;
[nb_ligne_1,nb_col_2]=size(max);


diff = diff(max(1,:));
mean_diff = mean(diff);
var_diff = var(diff);
k = 1;
for i=1:nb_col_2-1
    if((mean_diff-mean_diff/2)<=(max(1,i+1)-max(1,i)) & (max(1,i+1)-max(1,i))<=(mean_diff+mean_diff/2))
        corrected_max(1,k)=max(1,i);
        corrected_max(2,k)=max(2,i);
        k=k+1;
    end 
end
max = corrected_max;
[nb_ligne_1,nb_col_2]=size(max);



%% Recherche des points Q 

%1er minimum avant la R


for indice_4=1:nb_col_2
    loc_min_local = max(1,indice_4);
    while signal_de_base_retarde(1,loc_min_local)>=signal_de_base_retarde(1,loc_min_local-1)
        loc_min_local=loc_min_local-1;
    end
    point_Q(1,indice_4)=loc_min_local;
    point_Q(2,indice_4)=signal_de_base_retarde(1,loc_min_local);
    
end

        
%% Recherche des points S

%1er minimum apr�s la R


for indice_5=1:nb_col_2
    loc_min_local = max(1,indice_5);
    while signal_de_base_retarde(1,loc_min_local) >= signal_de_base_retarde(1,loc_min_local+1)
        loc_min_local = loc_min_local+1;
    end
    point_S(1,indice_5) = loc_min_local;
    point_S(2,indice_5) = signal_de_base_retarde(1,loc_min_local);
    
end

%% Superposition des points Q,R,S

figure(1)
hold on
title('point Q,R,S')
plot(signal_de_base_retarde)
stem(max(1,:),max(2,:));
stem(point_Q(1,:),point_Q(2,:));
stem(point_S(1,:),point_S(2,:));
hold off


%% Recherche des points T

G1_num=[1 0 0 0 0 0 -1];
G1_den=[1];

G2_num=[1 0 0 0 0 0 0 0 -1];
G2_den=[1 -1];

after_G1=filter(G1_num,G1_den,ecg);
after_G2=filter(G2_num,G2_den,after_G1);

%% On adapte les indices des points Q,R,S au signal de base 

point_R = max;

[nb_ligne_1,nb_col_2]=size(point_R);
for k=1:nb_col_2
    point_R(1,k) = point_R(1,k)-((N-1)/2+2+21); 
end 

[nb_ligne_1,nb_col_2]=size(point_S);
for k=1:nb_col_2
    point_S(1,k) = point_S(1,k)-((N-1)/2+2+21); 
end 

[nb_ligne_1,nb_col_2]=size(point_Q);
for k=1:nb_col_2
    point_Q(1,k) = point_Q(1,k)-((N-1)/2+2+21); 
end 


%% calcul du retard du pour trouver T

%fvtool(G1_num,G1_den) %nous donne un retard de 3
%fvtool(G2_num,G2_den) %donne environ 3.5


retard_3=zeros(1,14);
retard_3(1,3+4+1)=1; % diffrentiateur + low pass filter
after_G2_retard=conv(after_G2,retard_3); 
signal_de_base_retarde_2 = conv(ecg,retard_3);

[nb_ligne_1,nb_col_2]=size(point_R);
for k=1:nb_col_2
   point_R(1,k) = point_R(1,k)+7;
end 

[nb_ligne_1,nb_col_2]=size(point_S);
for k=1:nb_col_2
   point_S(1,k) = point_S(1,k)+7;
end 

[nb_ligne_1,nb_col_2]=size(point_Q);
for k=1:nb_col_2
   point_Q(1,k) = point_Q(1,k)+7;
end 


figure(2)
hold on
plot(signal_de_base_retarde_2)
plot(after_G2_retard)
stem(point_R(1,:),point_R(2,:));
stem(point_S(1,:),point_S(2,:));
stem(point_Q(1,:),point_Q(2,:));
hold off

            
[nb_ligne_1,nb_col_2]=size(point_S);

for i=1:nb_col_2-1
    indice_7=point_S(1,i);
    zero=0;
    indice_zero=1;
    while point_R(1,i)+0.6*(point_S(1,i+1)-point_S(1,i)) > indice_7 
        if ((after_G2_retard(1,indice_7)<0 & after_G2_retard(1,indice_7+1)>0 )| (after_G2_retard(1,indice_7)>0 & after_G2_retard(1,indice_7+1)<0))
          zero(1,indice_zero)=indice_7;
          zero(2,indice_zero)=signal_de_base_retarde_2(indice_7);
          indice_zero=indice_zero+1;
        end
        indice_7=indice_7+1;
    end
%     [nb_ligne_zero nb_col_zero]=size(zero);
%     max_zero=0;
%     max_zero=zero(2,1);
%     indice_max_zero=0;
%     for k=1:nb_col_zero
%         if (zero(2,k)>max_zero)
%             max_zero=zero(2,k);
%             indice_max_zero=zero(1,k);
%         end
%     end
%     point_T(1,i)=indice_max_zero;
%     point_T(2,i)=max_zero;
end




