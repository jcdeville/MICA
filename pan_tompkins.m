clear all;
close all;

load('ecg_normal_1.mat')
figure(1)
plot(abs(fft(ecg)))

%% part 1-2 : QRS complex - Pan & tompkins method


%low pass filter
low_pass_num=[1 0 0 0 0 0 -2 0 0 0 0 0 1];
low_pass_den=[1 -2 1];

after_lpf=filter(low_pass_num,low_pass_den,ecg);
figure(3)
hold on
plot(ecg)
plot(after_lpf)
%high pass filter
high_pass_num=zeros(1,33);
high_pass_num(1,1)=-1;
high_pass_num(1,17)=32;
high_pass_num(1,18)=-32;
high_pass_num(1,33)=1;

high_pass_den=[1 -1];

after_hpf=filter(high_pass_num,high_pass_den,after_lpf);
plot(after_hpf)

legend('ecg','after_lpf','after_hpf');

hold off
figure(2)
plot(abs(fft(after_hpf)))

%filter derivative

deri=[1 2 0 -2 -1].*Fs/8;
deriv=filter(deri,1,after_hpf);
figure(4)
plot(deriv)
%fvtool(deriv);

% squared signal

ssq=abs(deriv).^2;
figure(5)
plot(ssq)

%moving windows integration

N=0.10*Fs+1;


[nb_ligne,nb_col]=size(ssq);
for n=N:nb_col
    sum_1=0;
    for i=0:N-1
    sum_1=sum_1+ssq(1,n-i);
    end
    swni(1,n)=1/N*sum_1;
end

% tresholding

%ajout de 1 à N pour avoir un retard par car le retard est de N-1/2
%delai introduit par la window = 10, car en fait c'est une convolution avec une porte (la
%somme)

figure(6)
hold on

retard=zeros(1,14);
retard(1,(N-1)/2+2+1)=1; % retard windows + derivative
signal_retarder=conv(after_hpf,retard);
plot(swni)

plot(signal_retarder*10^8)
hold off

threshold=mean(swni); %mean of the signal

for j=1:nb_col
    if (swni(1,j)<threshold)
        swni(1,j)=0;
    end
end
figure(8)
hold on

plot(swni)

plot(signal_retarder*10^8)
hold off


retard_2=zeros(1,50);
retard_2(1,(N-1)/2+2+1+21)=1;

signal_de_base_retarde=conv(ecg,retard_2);

indice=1;
indice_3=1;


while indice<nb_col
     max_local=0;
     indice_2=0;
    indice_max=0;
    if(swni(1,indice)~=0)
        while(swni(1,indice+indice_2)~=0)
            if(max_local<signal_de_base_retarde(1,indice+indice_2))
               max_local=signal_de_base_retarde(1,indice+indice_2);
               indice_max=indice+indice_2;
            end
            indice_2=indice_2+1;
        end
        if(indice_max~=0)
         point_R(1,indice_3)=indice_max;
         point_R(2,indice_3)=signal_de_base_retarde(1,indice_max);
         indice_3=indice_3+1;

        end
        indice=indice+indice_2;
    else 
        indice=indice+1;
    end
end

[nb_ligne_1,nb_col_2]=size(point_R);

%idee 

mean_max=mean(point_R(2,:));
var_max=var(point_R(2,:));
k=1;
for i=1:nb_col_2
    if(point_R(2,i)>=mean_max-5*sqrt(var_max) & point_R(2,i)<=mean_max+5*sqrt(var_max))
        corrected_point_R(1,k)=point_R(1,i);
        corrected_point_R(2,k)=point_R(2,i);
                k=k+1;

    end
end


%superposer les points R sur ecg après passe_bande
close all

% figure()
% hold on
% plot(signal_retarder)
% stem(corrected_max(1,:),corrected_max(2,:));
% hold off

%mean heartbeat

[nb_ligne_3,nb_col_3]=size(corrected_point_R);


for h=2:nb_col_3
    gap_beatrate(1,h-1)=(corrected_point_R(1,h)-corrected_point_R(1,h-1));
end

mean_gap=mean(gap_beatrate);

beatrate=mean_gap/Fs;  %temps d'un battement en seconde
bpm=60/beatrate;

%% calcul du retard du filtre passe bande

%fvtool(low_pass_num,low_pass_den) %nous donne un retard de 5
%fvtool(high_pass_num,high_pass_den) %donne environ 16

%% superposion de ecg et des points R


figure()
title('point R')
hold on
plot(signal_de_base_retarde)
stem(point_R(1,:),point_R(2,:));
hold off

figure()
title('point R corrigé')
hold on
plot(signal_de_base_retarde)
stem(corrected_point_R(1,:),corrected_point_R(2,:));
hold off


%% recherche des points Q 

%1er minimum avant la R


for indice_4=1:nb_col_3
    loc_min_local=corrected_point_R(1,indice_4);
    while signal_de_base_retarde(1,loc_min_local)>=signal_de_base_retarde(1,loc_min_local-1)
        loc_min_local=loc_min_local-1;
    end
    point_Q(1,indice_4)=loc_min_local;
    point_Q(2,indice_4)=signal_de_base_retarde(1,loc_min_local);
    
end

        
figure()
hold on
title('point Q')
plot(signal_de_base_retarde)
stem(point_Q(1,:),point_Q(2,:));
hold off


%% recherche des points S

%1er minimum après la R


for indice_5=1:nb_col_3
    loc_min_local=corrected_point_R(1,indice_5);
    while signal_de_base_retarde(1,loc_min_local)>=signal_de_base_retarde(1,loc_min_local+1)
        loc_min_local=loc_min_local+1;
    end
    point_S(1,indice_5)=loc_min_local;
    point_S(2,indice_5)=signal_de_base_retarde(1,loc_min_local);
    
end

        
figure()
hold on
title('point S')
plot(signal_de_base_retarde)
stem(point_S(1,:),point_S(2,:));
hold off


%% part 1-2 : T and P points

% filter G1
G1_num=[1 0 0 0 0 0 -1];
% fvtool(G1_num,1) retard de 3

% filter G2
G2_num=[1 0 0 0 0 0 0 0 -1];
G2_den=[1 -1];

 after_G1=filter(G1_num,1,ecg);
 after_G2=filter(G2_num,G2_den,after_G1);

%fvtool(G2_num,G2_den) %retard de 3,5

%%
delay_G1_G2(1,8)=1;
signal_retard_3=conv(ecg,delay_G1_G2);

%% same delay as ecg

for i=1:nb_col_3
    point_S(1,i)=point_S(1,i)-((N-1)/2+2+21);
    point_R(1,i)=corrected_point_R(1,i)-((N-1)/2+2+21);
    point_R(2,i)=corrected_point_R(2,i);
    point_Q(1,i)=point_Q(1,i)-((N-1)/2+2+21);
end

%% same delay as ecg after G1 and G2 filter
for i=1:nb_col_3
    point_S_delay_G1_G2(1,i)=point_S(1,i)+7;
    point_S_delay_G1_G2(2,i)=point_S(2,i);
    point_R_delay_G1_G2(1,i)=point_R(1,i)+7;
    point_R_delay_G1_G2(2,i)=point_R(2,i);
    point_Q_delay_G1_G2(1,i)=point_Q(1,i)+7;
    point_Q_delay_G1_G2(2,i)=point_Q(2,i);

end

%% display ecg and the ecg after G1 and G2 filters 
figure()
hold on 
plot(signal_retard_3)
stem(point_S_delay_G1_G2(1,:),point_S_delay_G1_G2(2,:))
plot(after_G2*10^-1)
plot(0)
hold off

%% research of the T points

for i=1:nb_col_3-1
    j=point_S_delay_G1_G2(1,i);
    zero=0;
    indice_zero=1;
    while j< point_R_delay_G1_G2(1,i)+0.6*(point_S_delay_G1_G2(1,i+1)-point_S_delay_G1_G2(1,i))
        if ((after_G2(1,j)<0 & after_G2(1,j+1)>0 )| (after_G2(1,j)>0 & after_G2(1,j+1)<0))
        zero(1,indice_zero)=j;
        zero(2,indice_zero)=signal_retard_3(j);
        indice_zero=indice_zero+1;
        end
        j=j+1;
    end
         [nb_ligne_zero nb_col_zero]=size(zero);
         max_zero=0;
         max_zero=zero(2,1);
         indice_max_zero=0;
        for k=1:nb_col_zero
            if (zero(2,k)>max_zero)
                max_zero=zero(2,k);
                indice_max_zero=zero(1,k);
            end
    end
    point_T_delay_G1_G2(1,i)=indice_max_zero;
    point_T_delay_G1_G2(2,i)=max_zero;
    
end


%% research of the p points


for i=1:nb_col_3-1
    j=point_T_delay_G1_G2(1,i)+1;
    zero=0;
    indice_zero=1;
    while j<point_Q_delay_G1_G2(1,i+1)
        if ((after_G2(1,j)<0 & after_G2(1,j+1)>0 )| (after_G2(1,j)>0 & after_G2(1,j+1)<0))
        zero(1,indice_zero)=j;
        zero(2,indice_zero)=signal_retard_3(j);
        indice_zero=indice_zero+1;
        end
        j=j+1;
    end
         [nb_ligne_zero nb_col_zero]=size(zero);
         max_zero=0;
         max_zero=zero(2,1);
         indice_max_zero=zero(1,1);
        for k=1:nb_col_zero
            if (zero(2,k)>max_zero)
                max_zero=zero(2,k);
                indice_max_zero=zero(1,k);
            end
        end
    point_P_delay_G1_G2(1,i)=indice_max_zero;
    point_P_delay_G1_G2(2,i)=max_zero;   
end 
    


%% same delay as ecg
for i=1:nb_col_3-1
    point_T(1,i)=point_T_delay_G1_G2(1,i)-7;
    point_T(2,i)=point_T_delay_G1_G2(2,i);
    point_P(1,i)=point_P_delay_G1_G2(1,i)-7;
    point_P(2,i)=point_P_delay_G1_G2(2,i);
end
 
%% dipslay PQRST points

figure()
hold on
plot(ecg)
plot(point_P(1,:),point_P(2,:), '+','Color','red'); %text(point_P(1,:),point_P(2,:),' P ','Color','red','FontSize',14);
plot(point_Q(1,:),point_Q(2,:), '+','Color','red'); %text(point_Q(1,:),point_Q(2,:),' Q ','Color','red','FontSize',14);
plot(point_R(1,:),point_R(2,:), '+','Color','red'); %text(point_R(1,:),point_R(2,:),' R ','Color','red','FontSize',14);
plot(point_S(1,:),point_S(2,:), '+','Color','red'); %text(point_S(1,:),point_S(2,:),' S ','Color','red','FontSize',14);
plot(point_T(1,:),point_T(2,:), '+','Color','red'); %text(point_T(1,:),point_T(2,:),' T ','Color','red','FontSize',14);
hold off


