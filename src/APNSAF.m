clc
clear all
close all
sigmav=0.0988;
% Input: 
load coefficient_4_64;        % load filter bank parameters: M, hk, fk
Nr     = 500;
N      = 500;
B      = [0.32,-0.3,0.5,0.2];
AA      = 1; 
gamma  = 1e-2;
L      = M;                % interpolation/decimation factor
a      = 0.01;
Nw=6;
mu=0.4;

for ll=1:Nr
s=(2*randi([0,1],2000,1))-1+1j*(2*randi([0,1],2000,1)-1);

v=(sqrt(sigmav/2)*randn(2000,1))+(1j*(sqrt(sigmav/2)*randn(2000,1)));
for j=4:2000
X(j)=0.5*s(j)+1.2*s(j-1)+1.5*s(j-2)-s(j-3);
end
X(1)=0.5*s(1);
X(2)=0.5*s(2)+1.2*s(1);
X(3)=0.5*s(3)+1.2*s(2)+1.5*s(1);
d=filter(B,AA,X.');  
d=d+v;
desire1=zeros(4,1);
desire2=zeros(4,1);
desire3=zeros(4,1);
desire4=zeros(4,1);

for k=1:M
        uaux=filter(hk(k,:),1,d); %%d Hfilter output%%
	dsb(k,:)=uaux(find(mod((1:length(uaux))-1,L)==0)); % desired signal split in subbands
	uaux=filter(hk(k,:),1,X.');
	usb(k,:)=uaux(find(mod((1:length(uaux))-1,L)==0)); % input signal split in subbands
    end;
   
    for nn=10:500
   
    U1(:,nn)=usb(1,nn:-1:nn-Nw+1);
     U2(:,nn)=usb(2,nn:-1:nn-Nw+1);
      U3(:,nn)=usb(3,nn:-1:nn-Nw+1);
       U4(:,nn)=usb(4,nn:-1:nn-Nw+1);
    
    end

W1=zeros(6,500);
W2=zeros(6,500);
W3=zeros(6,500);
W4=zeros(6,500);

    for nn=100:500
for i=1:4
    A1(:,i)=U1(:,nn-i+1);
 A2(:,i)=U2(:,nn-i+1);
 A3(:,i)=U3(:,nn-i+1);
 A4(:,i)=U4(:,nn-i+1);
 
 desire1(i,1)=conj(dsb(1,nn-i+1));
desire2(i,1)=conj(dsb(2,nn-i+1));
desire3(i,1)=conj(dsb(3,nn-i+1));
desire4(i,1)=conj(dsb(4,nn-i+1));
  
end

W1(:,nn+1)=W1(:,nn)+mu*A1*(inv(A1'*A1))*(desire1(:,:)- A1'*W1(:,nn));
W2(:,nn+1)=W2(:,nn)+mu*A2*(inv(A2'*A2))*(desire2(:,:)- A2'*W2(:,nn));
W3(:,nn+1)=W3(:,nn)+mu*A3*(inv(A3'*A3))*(desire3(:,:)- A3'*W3(:,nn));
W4(:,nn+1)=W4(:,nn)+mu*A4*(inv(A4'*A4))*(desire4(:,:)- A4'*W4(:,nn));

e1f(:,nn)=(norm(dsb(1,nn)-W1(:,nn+1)'*A1(:,1))).^2;
e2f(:,nn)=(norm(dsb(2,nn)-W2(:,nn+1)'*A2(:,1))).^2;
e3f(:,nn)=(norm(dsb(3,nn)-W3(:,nn+1)'*A3(:,1))).^2;
e4f(:,nn)=(norm(dsb(4,nn)-W4(:,nn+1)'*A4(:,1))).^2;

e1(:,nn)=dsb(1,nn)-W1(:,nn+1)'*A1(:,1);
e2(:,nn)=dsb(2,nn)-W2(:,nn+1)'*A2(:,1);
e3(:,nn)=dsb(3,nn)-W3(:,nn+1)'*A3(:,1);
e4(:,nn)=dsb(4,nn)-W4(:,nn+1)'*A4(:,1);

e1_o=filter(fk(1,:),1,e1);
e2_o=filter(fk(2,:),1,e2);
e3_o=filter(fk(3,:),1,e3);
e4_o=filter(fk(4,:),1,e4);

efinal=e1_o+e2_o+e3_o+e4_o;
efinal2=abs((efinal)).^2;

  end
    
qq1(ll,:)=e1f;
qq2(ll,:)=e2f;
qq3(ll,:)=e3f;
qq4(ll,:)=e4f;
qqfin(ll,:)=efinal2;
end

    J1=(1./500)*(sum(qq1));
    J2=(1./500)*(sum(qq2));
    J3=(1./500)*(sum(qq3));
    J4=(1./500)*(sum(qq4));
    Jfin=(1./500)*(1./sum(qqfin));
    figure
    %plot(10*log10(J1))
    semilogy(J1)
    title('Learning Curve of error from subband 1st');
    xlabel('n');
    ylabel('MSE')
    grid on
    hold on
    figure
   % plot(10*log10(J2))
    semilogy(J2)
    axis([50,500,0.1,0.3])
    title('Learning Curve of error from subband 2nd');
    xlabel('n');
    ylabel('MSE')
    grid on
    hold on
    figure
   % plot(10*log10(J3))
   semilogy(J3)
   axis([50,500,0.1,1])
    title('Learning Curve of error from subband 3rd');
    xlabel('n');
    ylabel('MSE')
    grid on
    hold on
    figure
    semilogy(J4)
    title('Learning Curve of error from subband 4th');
    xlabel('n');
    ylabel('MSE')
    grid on
    hold on
    figure
    semilogy(Jfin)
    title('Learning Curve of total error');
    xlabel('n');
    ylabel('MSE')
    grid on
    hold on
    
