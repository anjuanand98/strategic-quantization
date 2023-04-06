function main
clc;
close all;
clear all;
%%
% X in range [a,b]
a=-5;
b=5;

%%
% discretizing theta, theta in [at,bt], mean mut, variance sigma_thsq
at=-5;
bt=5;
% thval=linspace(at+(bt-at)/(2*nt),bt-(bt-at)/(2*nt),nt); 
mut=0;
sigma_thsq=1;
thval1=linspace(at,mut-2*sigma_thsq,1);
thval2=linspace(mut-2*sigma_thsq,mut-sigma_thsq,2);
thval3=linspace(mut-sigma_thsq,mut+sigma_thsq,3);
thval4=linspace(mut+sigma_thsq,mut+2*sigma_thsq,2);
thval5=linspace(mut+2*sigma_thsq,bt,1);
thval=[thval1(2:end) thval2(2:end) thval3(2:end) thval4(2:end) thval5(2:end-1)];
thval=[thval1 thval2(2:end) thval3(2:end) thval4(2:end) thval5(2:end-1)];
nt=length(thval);
% pdf of theta
pth=zeros(1,length(thval));
f12=@(tv) ((1/sqrt(2*pi*sigma_thsq))*exp(-(tv-mut).^2/(2*sigma_thsq)));
sct=integral(f12,at,bt,'ArrayValued',true);
pth(1)=integral(f12,at,thval(1)+(thval(2)-thval(1))/2,'ArrayValued',true)/sct;
for i=2:length(thval)-1
    pth(i)=integral(f12,thval(i)-(thval(i)-thval(i-1))/2,thval(i)+(thval(i+1)-thval(i))/2,'ArrayValued',true)/sct;
end
pth(length(thval))=integral(f12,thval(end)-(thval(end)-thval(end-1))/2,bt,'ArrayValued',true)/sct;

%%
% parameters
rn=5; % # of random initializations
Mval=[2]; % # of levels of quantization
rhoval=[0 0.2  0.5 0.9]; % correlation
p_b=0; % bit error rate

%%
% initializing parameters to store values
endist=zeros(length(Mval),length(rhoval)); % encoder distortion for all M and correlation
dedist=zeros(length(Mval),length(rhoval)); % decoder distortion for all M and correlation
xmall=zeros(length(thval),Mval(end)+1,length(Mval),length(rhoval)); % quantizers for all M and correlation
exitm=zeros(length(Mval),length(rhoval)); % exitflag for all M and correlation
xinitm=zeros(length(thval),Mval(end)+1,length(Mval),length(rhoval)); % initialization used for all M and correlation
xrnall=zeros(length(thval),Mval(end)+1,length(Mval),rn,length(rhoval)); % quantizers for all M and correlation for all iterations
xrninitall=zeros(length(thval),Mval(end)+1,length(Mval),rn,length(rhoval)); % quantizer initializations for all M and correlation for all iterations
enall=zeros(length(Mval),rn,length(rhoval)); % encoder distortions for all M and correlation, all intializations

%%
% computing constraint matrices for all M
lt=length(thval);
A_all=zeros(lt*(Mval(end)-2),lt*(Mval(end)-1),length(Mval));
for M=Mval
    A=[];
    b1=[];
    if M>2
        A=zeros((M-1-1)*length(thval),(M-1)*length(thval));
        A1=[1 -1 zeros(1,M-1-2+(M-1)*(length(thval)-1))];
        i=0;
        for j=1:length(thval)
        A(i+1,:)=A1;
        i1=i+1;
        for i=(j-1)*(M-1-1)+2:(j-1)*(M-1-1)+M-2
            A(i,:)=circshift(A(i-1,:),1);
        end
        if length(i)==0
            i=i1;
        end
        A1=circshift(A(i,:),2);
        end
        
    end
    A;
    b1;
    A_all(1:lt*(M-2),1:lt*(M-1),find(M==Mval))=A;
end

%%
% deterministic initializations
rt=3; % # of deterministic initializations
x0init=zeros(length(thval),Mval(end)+1,length(Mval),rt); % deterministic initialization values

for m=1:length(Mval)
    if Mval(m)==2
        for j=1:length(thval)
        x0init(j,1:Mval(m)+1,m,1)=[a a+(b-a)*j/(10*length(thval)) b];
        x0init(j,1:Mval(m)+1,m,2)=[a b-(b-a)/(2*j) b];
        end
        x0init(:,1:Mval(m)+1,m,3)=repmat(linspace(a,b,Mval(m)+1),length(thval),1);
    else  
        for j=1:length(thval)
            xt1=linspace(a,a+(b-a)*j/(10*length(thval)),Mval(m));
            x0init(j,1:Mval(m)+1,m,1)=[a xt1(2:end) b];
            xt2=linspace(b-(b-a)/(2*j),b,Mval(m));
            x0init(j,1:Mval(m)+1,m,2)=[a xt2(1:end-1) b];   
        end
            xt3=linspace(a,b,Mval(m)+1);
            x0init(:,1:Mval(m)+1,m,3)=repmat(xt3,length(thval),1);
    end
end

%%
% main
for rhoind=[1:length(rhoval)] % loop over correlation
rho=rhoval(rhoind); % correlation value rho
mux=0; % mean of source X
sigma_xsq=1; % variance of source X

mux_corr=mux+rho*(sigma_xsq/sigma_thsq)^(1/2)*(thval(:)-mut); % mean of X conditional on theta 
sigma_xsq_corr=(1-rho^2)*sigma_xsq; % variance of X conditional on theta 
f1=@(xv,i) ((1/sqrt(2*pi*sigma_xsq_corr))*exp(-(xv-mux_corr(i)).^2/(2*sigma_xsq_corr)))*pth(i); % pdf of X conditional on theta

for M=Mval % loop over M
p_err=1-(1-p_b)^(log(M)/log(2)); % symbol error
c_1=p_err/(M-1); % 
c_2=1-M*c_1;
loopind=find(M==Mval); % loop index for M
% constraints
A=[];
b1=[];
if M>2
    A=A_all(1:lt*(M-2),1:lt*(M-1),(find(M==Mval))); % from the computation before
    b1=zeros(size(A,1),1); 
end
lb1=[a*ones(M-1,1)];
lb=repmat(lb1,length(thval),1);
ub1=[b*ones(M-1,1)];
ub=repmat(ub1,length(thval),1);

endist_rn=zeros(rn,1); % encoder distortions for each initialization
dedist_rn=zeros(rn,1); % decoder distortions for each initialization
xrn=zeros(length(thval),M+1,rn); % quantizer for given M for each initialization
xinitrn=zeros(length(thval),M+1,rn); % quantizer initializations for given M for each initialization
exitrn=zeros(rn,1); % exit flags for given M for each initialization
tic
for r=1:rn % loop over initializations
disp(strcat('M=',num2str(M),', rho=',num2str(rho),', iteration=',num2str(r),', bit error rate=',num2str(p_b))); % display all parameters 
if r<=rt % deterministic initializations
    x0=x0init(:,1:M+1,loopind,r);
else
    x0=[a*ones(length(thval),1) rand(length(thval),M-1)*(b-a)+a b*ones(length(thval),1)]; % random initializations
end
x1=x0;
x1=sort(x1')';
x0=sort(x0')'
x0=x0(:,2:end-1); % optimizing only decision values that are not boundaries
x0=x0';
x0=x0(:);
fun=@(x)f22fn(x,thval,f1,a,b,p_err,c_1,c_2,mux); % objective function
% options = optimoptions('fmincon','MaxFunctionEvaluations',90000000,'MaxIterations',90000000,'Display','iter','PlotFcn',{@optimplotx,@optimplotfval,@optimplotfirstorderopt});
options = optimoptions('fmincon','MaxFunctionEvaluations',90000000,'MaxIterations',90000000);

[x,fval,exitflag,output,lambda,grad,hessian]=fmincon(fun,x0,A,b1,[],[],lb,ub,[],options); % gradient descent


x=[a*ones(length(thval),1) reshape(x,M-1,length(thval))' b*ones(length(thval),1)]; % gradient descent output
exitrn(r)=exitflag % exit flag 
xm=x % quantizer 
ym=reconstruction(xm,thval,f1,p_err,c_1,c_2,mux) % reconstruction levels
[dist_enc]=encoderdistortion(xm,ym,f1,thval,p_err) % encoder distortion
[dist_dec]=decoderdistortion(xm,ym,f1,thval,p_err) % decoder distortion
endist_rn(r)=dist_enc; % encoder distortions for each initialization
dedist_rn(r)=dist_dec; % decoder distortions for each initialization
xrn(:,:,r)=xm; % quantizer for given M for each initialization
% xrnall(:,1:M+1,loopind,r,rhoind)=xm; % quantizers for all M and correlation for all iterations
xinitrn(:,:,r)=x1; % quantizer initializations for given M for each initialization
% xrninitall(:,1:M+1,loopind,r,rhoind)=x1; % quantizer initializations for all M and correlation for all iterations
enall(loopind,r,rhoind)=dist_enc; % encoder distortions for all M and correlation, all intializations
end
toc
for r=1:rn
    xrnall(:,1:M+1,loopind,r,rhoind)=xrn(:,:,r);
    xrninitall(:,1:M+1,loopind,r,rhoind)=xinitrn(:,:,r);
end

indl=find(exitrn==1); % finding valid runs of gradient descent using exit flag values
[in1,in2]=min(endist_rn(indl)); % finding minimum encoder distortion within the valid runs 
disp(strcat('Results for M=',num2str(M),', rho=',num2str(rho),', bit error rate=',num2str(p_b))); % display result for (M,rho) pair
disp('encoder and decoder distortions for all initializations:')
endist_rn(indl)
dedist_rn(indl)
disp('quantizer:')
x_opt=xrn(:,:,indl(in2(1))) % optimum quantizer for current (M,rho)
y_opt=reconstruction(x_opt,thval,f1,p_err,c_1,c_2,mux)
disp('encoder distortion:')
e_opt=endist_rn(indl(in2(1))) % optimum encoder distortion for current (M,rho)
disp('decoder distortion:')
d_opt=dedist_rn(indl(in2(1))) % optimum decoder distortion for current (M,rho)

endist(loopind,rhoind)=endist_rn(indl(in2(1))); % encoder distortion for all M and correlation
dedist(loopind,rhoind)=dedist_rn(indl(in2(1))); % decoder distortion for all M and correlation
xmall(:,1:M+1,loopind,rhoind)=xrn(:,:,indl(in2(1))); % quantizers for all M and correlation
exitm(loopind,rhoind)=exitrn(indl(in2(1))); % exitflag for all M and correlation
xinitm(:,1:M+1,loopind,rhoind)=xrninitall(:,1:M+1,loopind,indl(in2(1))); % initialization used for all M and correlation
save(strcat('fmincon1_data_nt',num2str(nt),'noisy',num2str(p_b),'rho',num2str(rho),'M',num2str(M),'.mat'),'rho','p_b','x_opt','y_opt','e_opt','d_opt','thval','M','xmall','endist','dedist','exitm','xinitm','xrnall','xrninitall','enall');
end

% plot figure: encoder and decoder distortion for a given bit error p_b and correlation rho

% % f=figure;
% % plot(log(Mval)./log(2),endist,'-o')
% % hold on;
% % plot(log(Mval)./log(2),dedist,'-*')
% % hold off;
% % grid on;
% % xlabel('rate (in bits)')
% % ylabel('distortion')
% % legend('encoder distortion','decoder distortion')
% % saveas(f,strcat('error',num2str(p_b),'rho',num2str(rho),'theta5randinit.fig'))
% % saveas(f,strcat('error',num2str(p_b),'rho',num2str(rho),'theta5randinit.png'))

% plot figure: quantizer for each M for given bit error p_b and correlation rho

% % for i=1:length(Mval)
% % f=figure;
% % for j=1:length(thval)
% %     plot(xmall(j,1:Mval(i)+1,i),thval(j)*ones(1,Mval(i)+1),'-o')
% %     hold on;
% % end
% % hold off;
% % grid on;
% % ylim([at bt])
% saveas(f,strcat('error',num2str(p_b),'rho',num2str(rho),'quantizer',num2str(Mval(i)),'_',num2str(nt),'.png'))
% saveas(f,strcat('error',num2str(p_b),'rho',num2str(rho),'quantizer',num2str(Mval(i)),'_',num2str(nt),'.fig'))
% % end

end
save(strcat('fmincon1_data_nt',num2str(nt),'noisy',num2str(p_b),'rho.mat'),'thval','rhoval','Mval','p_b','xmall','endist','dedist','exitm','xinitm','xrnall','xrninitall','enall');

function [ym]=reconstruction(xthetam,thval,f1,p_err,c_1,c_2,mux)
M=size(xthetam,2)-1;
ym=zeros(1,M);
for i=1:M
    num=0;
    den=0;
    for j=1:length(thval)
        f1temp= @(xv) f1(xv,j);
        f2=@(xv) xv*f1temp(xv);
        num=num+integral(f2,xthetam(j,i),xthetam(j,i+1),'ArrayValued',true);
        den=den+integral(f1temp,xthetam(j,i),xthetam(j,i+1),'ArrayValued',true);
    end
    if den~=0
        ym(i)=(c_1*mux+c_2*num)/(c_1+c_2*den);
    else
        ym(i)=(1/size(xthetam,1))*sum(xthetam(:,i));
    end
end

function [dist_dec]=decoderdistortion(xthetam,ym,f1,thval,p_err)
M=size(xthetam,2)-1;
c1=p_err/(M-1);
c2=1-M*c1;
dist_dec=0;
for j=1:M%yj
    for i=1:M
        for k=1:length(thval)
            f1temp= @(xv) f1(xv,k);
            f5=@(xv) (xv-ym(j))^2*f1temp(xv);
            dist_dec=dist_dec+c1*integral(f5,xthetam(k,i),xthetam(k,i+1),'ArrayValued',true);
        end
    end
end
for i=1:M
    for k=1:length(thval)
        f1temp= @(xv) f1(xv,k);
        f5=@(xv) (xv-ym(i))^2*f1temp(xv);
        dist_dec=dist_dec+c2*integral(f5,xthetam(k,i),xthetam(k,i+1),'ArrayValued',true);
    end
end

function [dist_enc]=encoderdistortion(xthetam,ym,f1,thval,p_err)
M=size(xthetam,2)-1;
c1=p_err/(M-1);
c2=1-M*c1;
dist_enc=0;
for j=1:M%yj
    for i=1:M
        for k=1:length(thval)
            f1temp= @(xv) f1(xv,k);
            f5=@(xv) (xv+thval(k)-ym(j))^2*f1temp(xv);
            dist_enc=dist_enc+c1*integral(f5,xthetam(k,i),xthetam(k,i+1),'ArrayValued',true);
        end
    end
end
for i=1:M
    for k=1:length(thval)
        f1temp= @(xv) f1(xv,k);
        f5=@(xv) (xv+thval(k)-ym(i))^2*f1temp(xv);
        dist_enc=dist_enc+c2*integral(f5,xthetam(k,i),xthetam(k,i+1),'ArrayValued',true);
    end
end

function [f22] = f22fn(x,thval,f1,a,b,p_err,c_1,c_2,mux)

M=length(x)/length(thval)+1;
x=[a*ones(length(thval),1) reshape(x,M-1,length(thval))' b*ones(length(thval),1)];
[ym]=reconstruction(x,thval,f1,p_err,c_1,c_2,mux);
x=x';
x=x(:);

f22=0;
for i=1:M
    for t=1:length(thval)
            f22=f22+c_2*integral(@(xv)(xv+thval(t)-ym(i))^2*f1(xv,t),x((t-1)*(M+1)+i),x((t-1)*(M+1)+i+1),'ArrayValued',true);
        for j=1:M
            f22=f22+c_1*integral(@(xv)(xv+thval(t)-ym(j))^2*f1(xv,t),x((t-1)*(M+1)+i),x((t-1)*(M+1)+i+1),'ArrayValued',true);
        end
    end
end

