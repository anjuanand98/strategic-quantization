% figures: noisy gaussian, (x^3-y)^2
clc;
close all;
clear all;
pbval=[0 0.05 0.2 0.49];
Mval=[2 4 8];
load(strcat('xcubed_noisy_fmincon1_data'));
stle=strings(1,length(pbval));
f=figure;
    plot(log(Mval)/log(2),reshape(endist(:,1),1,length(Mval)),'--o','LineWidth',2.0,'Markersize',15);
    stle(1)=strcat('p_{b}=',num2str(pbval(1)));
    hold on;
    plot(log(Mval)/log(2),reshape(endist(:,2),1,length(Mval)),'--^','LineWidth',2.0,'Markersize',15);
    stle(2)=strcat('p_{b}=',num2str(pbval(2)));
    hold on;
    plot(log(Mval)/log(2),reshape(endist(:,3),1,length(Mval)),'--d','LineWidth',2.0,'Markersize',15);
    stle(3)=strcat('p_{b}=',num2str(pbval(3)));
    hold on;
    plot(log(Mval)/log(2),reshape(endist(:,4),1,length(Mval)),'--p','LineWidth',2.0,'Markersize',15);
    stle(4)=strcat('p_{b}=',num2str(pbval(4)));
    hold off;
    grid on;
    lgd=legend(stle);
    lgd.FontSize=14;
    lgd.NumColumns=2;
    xlabel('rate (in bits)','FontSize',14)
    ylabel('encoder distortion','FontSize',14)
    saveas(f,strcat('encdist_xcubed_noisy.png'))
    saveas(f,strcat('encdist_xcubed_noisy.fig'))
    f=figure;
    plot(log(Mval)/log(2),reshape(dedist(:,1),1,length(Mval)),'--o','LineWidth',2.0,'Markersize',15);
    stle(1)=strcat('p_{b}=',num2str(pbval(1)));
    hold on;
    plot(log(Mval)/log(2),reshape(dedist(:,2),1,length(Mval)),'--^','LineWidth',2.0,'Markersize',15);
    stle(2)=strcat('p_{b}=',num2str(pbval(2)));
    hold on;
    plot(log(Mval)/log(2),reshape(dedist(:,3),1,length(Mval)),'--d','LineWidth',2.0,'Markersize',15);
    stle(3)=strcat('p_{b}=',num2str(pbval(3)));
    hold on;
    plot(log(Mval)/log(2),reshape(dedist(:,4),1,length(Mval)),'--p','LineWidth',2.0,'Markersize',15);
    stle(4)=strcat('p_{b}=',num2str(pbval(4)));
    hold off;
    grid on;
    lgd=legend(stle);
    lgd.FontSize=14;
    lgd.NumColumns=2;
    xlabel('rate (in bits)','FontSize',14)
    ylabel('decoder distortion','FontSize',14)
    saveas(f,strcat('decdist_xcubed_noisy.png'))
    saveas(f,strcat('decdist_xcubed_noisy.fig'))