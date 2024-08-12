clc;clear all; 
close all
delta=15;
N=4;
MM=6;
mu=0.4;
n=500;
W=zeros(MM,500);
sigmav=0.0988;
B      = [0.32,-0.3,0.5,0.2];
AA      = 1; 
load coefficient_4_64;
for jj=1:500
s=(2*randi([0,1],n,1))-1+1j*(2*randi([0,1],n,1)-1);

v=(sqrt(sigmav/2)*randn(n,1))+(1j*(sqrt(sigmav/2)*randn(n,1)));
for j=4:n
X(j)=0.5*s(j)+1.2*s(j-1)+1.5*s(j-2)-s(j-3);
end
X(1)=0.5*s(1);
X(2)=0.5*s(2)+1.2*s(1);
X(3)=0.5*s(3)+1.2*s(2)+1.5*s(1);
U=X.'+v;
desire=filter(B,AA,X.');  
desire=desire+v;
d=zeros(N,1);

for nn=10:500
   
    u(:,nn)=U(nn:-1:nn-MM+1);
    
end
for nn=100:500
for i=1:N
    A(:,i)=u(:,nn-i+1);
    d(i,1)=conj(desire(nn-i+1,1));
  
end
 %d1=s(nn-15);
W(:,nn+1)=W(:,nn)+mu*A*(inv(A'*A))*(d-A'*W(:,nn));
e(:,nn)=(norm(desire(nn)-W(:,nn+1)'*A(:,1))).^2;
% e1=e(nn);
end
qq(jj,:)=(e);
end
J=(1./500)*(sum(qq))
semilogy(J)
title('Learning Curve of APA for N=4');
xlabel('n');
ylabel('MSE')
grid on
hold on

