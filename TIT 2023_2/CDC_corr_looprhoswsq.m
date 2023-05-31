function main
clc;
close all;
clear all;

%%
rhoval=[-0.99 -0.2];
sigma_wsqval=[0.01 0.2 0.6];
%%
au=-5;
bu=5;
aw=-5;
bw=5;
% rho=-0.2;
%%
% discretizing theta, theta in [at,bt], mean mut, variance sigma_thsq
at=-5;
bt=5;
% X in range [ax,bx,au,bu,]
ax=au+aw;
bx=bu+bw;
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
Mval=[ 2 4 ];
endistM=zeros(length(Mval),1);
dedistM=zeros(length(Mval),1);
for sigma_wsq=sigma_wsqval
for rho=rhoval
for M=Mval
%%
% computing constraint matrices 

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
        b1=zeros(size(A,1),1); 
    end

% forming a grid
samp=linspace(-3+0.001,3-0.001,15);
samp=[-8 -6 -4 samp 4 6 8];
x0init=nchoosek(samp(2:end-1),M-1);
x0init=[ax*ones(size(x0init,1),1) x0init bx*ones(size(x0init,1),1)];
rn=20;
xrn1=randi(size(x0init,1),rn,length(thval));
x0r=zeros(rn,length(thval),M+1);
xmoptrn=zeros(rn,length(thval),M+1);
firstordoptrn=zeros(rn,1);
for r=1:rn

x0r(r,:,:)=x0init(xrn1(r,:)',:);
end
edistr=zeros(rn,1);
ddistr=zeros(rn,1);
xdistr=zeros(rn,length(thval),M+1);
exitr=zeros(rn,1);
%%
% main
muu=0; % mean of noiseless source U
sigma_usq=1; % variance of noiseless source U
muu_corr=muu+rho*(sigma_usq/sigma_thsq)^(1/2)*(thval(:)-mut); % mean of X conditional on theta 
sigma_usq_corr=(1-rho^2)*sigma_usq; % variance of X conditional on theta 
f1=@(uv,i) ((1/sqrt(2*pi*sigma_usq_corr))*exp(-(uv-muu_corr(i)).^2/(2*sigma_usq_corr)))*pth(i); % pdf of X conditional on theta


fu=@(uv) ((1/sqrt(2*pi*sigma_usq))*exp(-(uv-muu).^2/(2*sigma_usq))); % pdf of U 
scalu=integral(@(uv) fu(uv),au,bu);
fu=@(uv) fu(uv)./scalu;
muw=0; % mean of independent additive noise
% sigma_wsq=0.1; % variance of independent additive noise

fw=@(wv) ((1/sqrt(2*pi*sigma_wsq))*exp(-(wv-muw).^2/(2*sigma_wsq))); % pdf of W 
scalw=integral(@(wv) fw(wv),aw,bw);
fw=@(wv) fw(wv)/scalw;
% noisy source x=u+w
mux=muu+muw;
sigma_xsq=sigma_usq+sigma_wsq;
fx=@(xv) ((1/sqrt(2*pi*sigma_xsq))*exp(-(xv-mux).^2/(2*sigma_xsq))); % pdf of X
scalx=integral(@(xv) fx(xv),ax,bx);
fx=@(xv) fx(xv)/scalx;
fux=@(uv,xv) fu(uv).*fw(xv-uv);
fuxt=@(uv,xv,tv) f1(uv,tv).*fw(xv-uv);
% parameters used for computing E(U|x)
m1=muu-sigma_usq* (sigma_usq+sigma_wsq)^(-1)*(muu+muw);
K1=(sigma_usq)*(sigma_usq+sigma_wsq)^(-1);
lb1=[ax*ones(M-1,1)];
lb=repmat(lb1,length(thval),1);
ub1=[bx*ones(M-1,1)];
ub=repmat(ub1,length(thval),1);

parfor r=1:rn
x0=reshape(x0r(r,:,:),length(thval),M+1);
% x0=[a*ones(length(thval),1) rand(length(thval),M-1)*(b-a)+a b*ones(length(thval),1)]; % random initializations
x1=x0;
x1=sort(x1')';
x0=sort(x0')';
x0=x0(:,2:end-1); % optimizing only decision values that are not boundaries
x0=x0';
x0=x0(:);

fun=@(x)f22fn(x,thval,ax,bx,au,bu,fux,m1,K1,pth,fx,fuxt); % objective function
% options = optimoptions('fmincon','MaxFunctionEvaluations',90000000,'MaxIterations',90000000,'Display','iter','PlotFcn',{@optimplotx,@optimplotfval,@optimplotfirstorderopt});
options = optimoptions('fmincon','MaxFunctionEvaluations',90000000,'MaxIterations',90000000,'Display','final-detailed');
% tic;
[x,fval,exitflag,output,lambda,grad,hessian]=fmincon(fun,x0,A,b1,[],[],lb,ub,[],options); % gradient descent
% toc;
firstordoptrn(r)=output.firstorderopt;
x11=[ax*ones(length(thval),1) reshape(x,M-1,length(thval))' bx*ones(length(thval),1)]; % gradient descent output
exitflag % exit flag 
xm=x11 % quantizer 
xmoptrn(r,:,:)=x11;
ym=reconstruction(x11,thval,fux,pth,ax,bx,au,bu,fuxt) % reconstruction levels
[encoder_dist] = f22fn(x,thval,ax,bx,au,bu,fux,m1,K1,pth,fx,fuxt)% PLUS SOMETHING
[decoder_dist]=decoderdistortion(x11,ym,fux,pth,ax,bx,au,bu,fuxt)
exitr(r)=exitflag;
edistr(r)=encoder_dist;
ddistr(r)=decoder_dist;
xdistr(r,:,:)=xm;
end

indl=find(exitr>=1 & firstordoptrn<10^-3);
[m11,m12]=min(edistr(indl));
endistM(find(M==Mval))=edistr(indl(m12(1)));
encoder_dist=endistM(find(M==Mval));
x11=reshape(xdistr(indl(m12),:,:),length(thval),M+1);
ym=reconstruction(x11,thval,fux,pth,ax,bx,au,bu,fuxt);
[decoder_dist]=decoderdistortion(x11,ym,fux,pth,ax,bx,au,bu,fuxt);
dedistM(find(M==Mval))=decoder_dist;
save(strcat('M',num2str(M),'wvar',num2str(sigma_wsq),'rho',num2str(rho),'data.mat'),'encoder_dist','decoder_dist','x11','ym','exitr','edistr','x0r','xmoptrn','firstordoptrn','rho','sigma_wsq');
end
end
end
save(strcat('wvar',num2str(sigma_wsq),'rho',num2str(rho),'data.mat'),'endistM','dedistM');

function [dist_dec]=decoderdistortion(xthetam,ym,fux,pth,ax,bx,au,bu,fuxt)
M=size(xthetam,2)-1;
dist_dec=0;
for i=1:M
    for k=1:length(pth)
        fux2=@(uv,xv) (uv-ym(i)).^2.*fuxt(uv,xv,k);
        dist_dec=dist_dec+integral2(fux2,au,bu,xthetam(k,i),xthetam(k,i+1));
    end
end

function [ym]=reconstruction(xthetam,thval,fux,pth,ax,bx,au,bu,fuxt)
M=size(xthetam,2)-1;
ym=zeros(1,M);
for i=1:M
    num=0;
    den=0;
    for j=1:length(thval)
        fux1= @(uv,xv) uv.*fuxt(uv,xv,j);
        num=num+integral2(fux1,au,bu,xthetam(j,i),xthetam(j,i+1));
        den=den+integral2(@(uv,xv) fuxt(uv,xv,j),au,bu,xthetam(j,i),xthetam(j,i+1));
    end
    ym(i)=num/den;
end

function [f22] = f22fn(x,thval,ax,bx,au,bu,fux,m1,K1,pth,fx,fuxt)
M=length(x)/length(thval)+1;
x=[ax*ones(length(thval),1) reshape(x,M-1,length(thval))' bx*ones(length(thval),1)];
[ym]=reconstruction(x,thval,fux,pth,ax,bx,au,bu,fuxt);
x=x';
x=x(:);

% f22=0;
% for i=1:M
%     for t=1:length(thval)
%             f22=f22+integral(@(xv)((m1+K1*xv)+thval(t)-ym(i))^2*fx(xv)*pth(t),x((t-1)*(M+1)+i),x((t-1)*(M+1)+i+1),'ArrayValued',true);
%     end
% end
f22=0;
for i=1:M
    for t=1:length(thval)
            f22=f22+integral2(@(uv,xv)(uv+thval(t)-ym(i)).^2.*fuxt(uv,xv,t),au,bu,x((t-1)*(M+1)+i),x((t-1)*(M+1)+i+1));
    end
end

