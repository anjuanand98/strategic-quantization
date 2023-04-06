% figures: given rho, wrt pb
clc;
close all;
clear all;
rhoval=[0.99];
pbval=[0 0.05 0.2 0.49];
Mval=[2 4 8];
enc_dist=zeros(length(rhoval),length(pbval),length(Mval));
dec_dist=zeros(length(rhoval),length(pbval),length(Mval));

non_strdist=zeros(length(pbval),length(Mval));
for rho=rhoval
    rind=find(rho==rhoval);
for pb=pbval
    pbind=find(pb==pbval);
    for M=Mval
        Mind=find(M==Mval);
        load(strcat('fmincon1_data_nt5noisy',num2str(pb),'rho',num2str(rho),'M',num2str(M),'.mat'));
        enc_dist(rind,pbind,Mind)=e_opt;
        dec_dist(rind,pbind,Mind)=d_opt;
    end
end
end

for pb=pbval
    pbind=find(pb==pbval);
    for M=Mval
        Mind=find(M==Mval);
        load(strcat('nonstr_fmincon1_data_nt5noisy',num2str(pb),'rho',num2str(rho),'M',num2str(M),'.mat'));
        non_strdist(pbind,Mind)=d_opt;
    end
end


stle=strings(1,length(pbval));
for rho=rhoval
    f=figure;
    rind=find(rho==rhoval);
    plot(log(Mval)/log(2),reshape(dec_dist(rind,1,:),1,length(Mval)),'--o','LineWidth',2.0,'Markersize',15);
    stle(1)=strcat('p_{b}=',num2str(pbval(1)));
    hold on;
    plot(log(Mval)/log(2),reshape(dec_dist(rind,2,:),1,length(Mval)),'--^','LineWidth',2.0,'Markersize',15);
    stle(2)=strcat('p_{b}=',num2str(pbval(2)));
    hold on;
    plot(log(Mval)/log(2),reshape(dec_dist(rind,3,:),1,length(Mval)),'--d','LineWidth',2.0,'Markersize',15);
    stle(3)=strcat('p_{b}=',num2str(pbval(3)));
    hold on;
    plot(log(Mval)/log(2),reshape(dec_dist(rind,4,:),1,length(Mval)),'--p','LineWidth',2.0,'Markersize',15);
    stle(4)=strcat('p_{b}=',num2str(pbval(4)));
    hold on;

%     plot(log(Mval)/log(2),non_strdist(1,:))
%     hold on;
%     plot(log(Mval)/log(2),non_strdist(2,:))
%     hold on
%     plot(log(Mval)/log(2),non_strdist(3,:))
%     hold on
%     plot(log(Mval)/log(2),non_strdist(4,:))
    hold off;
    grid on;
    ax = gca;
    ax.FontSize = 18;
    lgd=legend(stle);
    lgd.FontSize=18;
    lgd.NumColumns=1;
    xlabel('rate (in bits)','FontSize',20)
    ylabel('strategic decoder distortion','FontSize',20)
%     saveas(f,strcat('encdist_rho',num2str(rho),'.png'))
%     saveas(f,strcat('encdist_rho',num2str(rho),'.fig'))
    f=figure;
    pbind=1;
    plot(log(Mval)/log(2),reshape(non_strdist(pbind,:),1,length(Mval)),'--o','LineWidth',2.0,'Markersize',15);
    stle(1)=strcat('p_{b}=',num2str(pbval(1)));
    hold on;
    pbind=2;
    plot(log(Mval)/log(2),reshape(non_strdist(pbind,:),1,length(Mval)),'--^','LineWidth',2.0,'Markersize',15);
    stle(2)=strcat('p_{b}=',num2str(pbval(2)));
    hold on;
    pbind=3;
    plot(log(Mval)/log(2),reshape(non_strdist(pbind,:),1,length(Mval)),'--d','LineWidth',2.0,'Markersize',15);
    stle(3)=strcat('p_{b}=',num2str(pbval(3)));
    hold on;
    pbind=4;
    plot(log(Mval)/log(2),reshape(non_strdist(pbind,:),1,length(Mval)),'--p','LineWidth',2.0,'Markersize',15);
    stle(4)=strcat('p_{b}=',num2str(pbval(4)));
    hold off;
    grid on;
    lgd=legend(stle);
    lgd.FontSize=18;
    lgd.NumColumns=1;
    ax = gca;
    ax.FontSize = 18;
    xlabel('rate (in bits)','FontSize',20)
    ylabel('non-strategic decoder distortion','FontSize',20)
    saveas(f,strcat('decdist_rho',num2str(rho),'.png'))
    saveas(f,strcat('decdist_rho',num2str(rho),'.fig'))
end