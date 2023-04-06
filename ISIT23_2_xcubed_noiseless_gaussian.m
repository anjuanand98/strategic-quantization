function main
clc;
close all;
clear all;

%%
% X in range [a,b] with probability distribution f1
a=-5;
b=5;
mux=0;
sigma_xsq=1;
f1=@(xv) ((1/sqrt(2*pi*sigma_xsq))*exp(-(xv-mux).^2/(2*sigma_xsq)));

%%
% parameters
rn=5;
Mval=[2 4 8];

%%
% initializing parameters to store values
endist=zeros(length(Mval),1); % encoder distortion for all M 
dedist=zeros(length(Mval),1); % decoder distortion for all M 
xmall=zeros(Mval(end)+1,length(Mval)); % quantizers for all M  
exitm=zeros(length(Mval),1); % exitflag for all M
xinitm=zeros(Mval(end)+1,length(Mval)); % initialization used for all M 
xrnall=zeros(Mval(end)+1,length(Mval),rn); % quantizers for all M for all iterations
xrninitall=zeros(Mval(end)+1,length(Mval),rn); % quantizer initializations for all M for all iterations
enall=zeros(length(Mval),rn); % encoder distortions for all M, all intializations

%%
% computing constraint matrices for all M
A_all=zeros(Mval(end)-2,Mval(end)-1,length(Mval));
for M=Mval
A=[];
b1=[];
if M>2
A= zeros(M-1-1,M-1);
A(1,:)=[1 -1 zeros(1,M-1-2)];
for i=2:M-1-1
    A(i,:)=circshift(A(i-1,:),1);
end
b1=zeros(size(A,1),1);
end
A_all(1:(M-2),1:(M-1),find(M==Mval))=A;
end

%% 
% deterministic initializations
rt=3;
x0init=zeros(Mval(end)+1,length(Mval),rt); % deterministic initialization values

for m=1:length(Mval)
    if Mval(m)==2
        x0init(1:Mval(m)+1,m,1)=[a; a+(b-a)/10; b];
        x0init(1:Mval(m)+1,m,2)=[a; b-(b-a)/2; b];
        x0init(1:Mval(m)+1,m,3)=[a; (a+b)/2; b];
    else
        xt1=linspace(a,a+(b-a)/10,Mval(m));
        x0init(1:Mval(m)+1,m,1)=[a; xt1(2:end)'; b];
        xt2=linspace(b-(b-a)/2,b,Mval(m));
        x0init(1:Mval(m)+1,m,2)=[a; xt1(1:end-1)'; b];
        xt3=linspace(a,b,Mval(m)+1);
        x0init(1:Mval(m)+1,m,3)=xt3';
    end
end

%%
% main
for M=Mval % loop over M
    loopind=find(M==Mval); % loop index for M
    % constraints
    A=[];
    b1=[];
    if M>2
        A=A_all(1:(M-2),1:(M-1),(find(M==Mval))); % from the computation before
        b1=zeros(size(A,1),1); 
    end
    lb=[a*ones(M-1,1)];
    ub=[b*ones(M-1,1)];
    endist_rn=zeros(rn,1); % encoder distortions for each initialization
    dedist_rn=zeros(rn,1); % decoder distortions for each initialization
    xrn=zeros(M+1,rn); % quantizer for given M for each initialization
    xinitrn=zeros(M+1,rn); % quantizer initializations for given M for each initialization
    exitrn=zeros(rn,1); % exit flags for given M for each initialization
    tic
    for r=1:rn % loop over initializations
        disp(strcat('M=',num2str(M),', iteration=',num2str(r))); % display all parameters 
        if r<=rt % deterministic initializations
            x0=x0init(1:M+1,loopind,r);
        else
            x0=[a; rand(M-1,1)*(b-a)+a; b]; % random initializations
        end
        x1=x0;
        x1=sort(x1')';
        x0=sort(x0')'
        x0=x0(2:end-1); % optimizing only decision values that are not boundaries
        fun=@(x)f22fn(x,f1,a,b); % objective function
        options = optimoptions('fmincon','MaxFunctionEvaluations',90000000,'MaxIterations',90000000);
        [x,fval,exitflag,output,lambda,grad,hessian]=fmincon(fun,x0,A,b1,[],[],lb,ub,[],options); % gradient descent
        x=[a; x; b]; % gradient descent output
        exitrn(r)=exitflag % exit flag 
        xm=x % quantizer 
        ym=reconstruction(xm,f1) % reconstruction levels
        [dist_enc]=encoderdistortion(xm,ym,f1) % encoder distortion
        [dist_dec]=decoderdistortion(xm,ym,f1) % decoder distortion
        endist_rn(r)=dist_enc; % encoder distortions for each initialization
        dedist_rn(r)=dist_dec; % decoder distortions for each initialization
        xrn(:,r)=xm; % quantizer for given M for each initialization
        % xrnall(:,1:M+1,loopind,r,rhoind)=xm; % quantizers for all M and correlation for all iterations
        xinitrn(:,r)=x1; % quantizer initializations for given M for each initialization
        % xrninitall(:,1:M+1,loopind,r,rhoind)=x1; % quantizer initializations for all M and correlation for all iterations
        enall(loopind,r)=dist_enc; % encoder distortions for all M and correlation, all intializations
    end
    toc
    for r=1:rn
    xrnall(1:M+1,loopind,r)=xrn(:,r);
    xrninitall(1:M+1,loopind,r)=xinitrn(:,r);
    end
    indl=find(exitrn==1); % finding valid runs of gradient descent using exit flag values
    [in1,in2]=min(endist_rn(indl)); % finding minimum encoder distortion within the valid runs 
    disp(strcat('Results for M=',num2str(M))); % display result for M
    disp('encoder and decoder distortions for all initializations:')
    endist_rn(indl)
    dedist_rn(indl)
    disp('quantizer:')
    x_opt=xrn(:,indl(in2(1))) % optimum quantizer for current M
    y_opt=reconstruction(x_opt,f1)
    disp('encoder distortion:')
    e_opt=endist_rn(indl(in2(1))) % optimum encoder distortion for current M
    disp('decoder distortion:')
    d_opt=dedist_rn(indl(in2(1))) % optimum decoder distortion for current M

    endist(loopind)=endist_rn(indl(in2(1))); % encoder distortion for all M 
    dedist(loopind)=dedist_rn(indl(in2(1))); % decoder distortion for all M 
    xmall(1:M+1,loopind)=xrn(:,indl(in2(1))); % quantizers for all M 
    exitm(loopind)=exitrn(indl(in2(1))); % exitflag for all M 
    xinitm(1:M+1,loopind)=xrninitall(1:M+1,loopind,indl(in2(1))); % initialization used for all M 
    save(strcat('xcubed_noiseless_gaussian_fmincon1_data','M',num2str(M),'.mat'),'x_opt','y_opt','e_opt','d_opt','M','xmall','endist','dedist','exitm','xinitm','xrnall','xrninitall','enall');
end
save(strcat('xcubed_noiseless_gaussian_fmincon1_data.mat'),'Mval','xmall','endist','dedist','exitm','xinitm','xrnall','xrninitall','enall');

function [ym]=reconstruction(xm,f1)
M=length(xm)-1;
f2=@(xv) xv*f1(xv);
ym=zeros(1,M);
for i=1:M
    if xm(i)~=xm(i+1)
    ym(i)=integral(f2,xm(i),xm(i+1),'ArrayValued',true)/integral(f1,xm(i),xm(i+1),'ArrayValued',true);
    else
        ym(i)=xm(i);
    end
end

function [dist_dec]=decoderdistortion(xm,ym,f1)
M=length(xm)-1;
dist_dec=0;
for i=1:M
    f5=@(xv) (xv-ym(i))^2*f1(xv);
    dist_dec=dist_dec+integral(f5,xm(i),xm(i+1),'ArrayValued',true);
end

function [dist_enc]=encoderdistortion(xm,ym,f1)
M=length(xm)-1;
dist_enc=0;
for i=1:M
    f4=@(xv) (xv^3-ym(i))^2*f1(xv);
    dist_enc=dist_enc+integral(f4,xm(i),xm(i+1),'ArrayValued',true);
end


function [f22] = f22fn(x,f1,a,b)

M=length(x)+1;
x=[a;x;b];
[ym]=reconstruction(x,f1);
    
f22=0;
for m=1:M
    f22=f22+ integral(@(xv)(xv^3-ym(m))^2*f1(xv),x(m),x(m+1),'ArrayValued',true);
end