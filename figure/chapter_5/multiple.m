%clear all;
close all;
N=16;
dt=0.1;
step=100;
d0=1;
d1=2.25;
K=1;
sigma=1;
beta=0.25;
theta=2;
avg=zeros(1,step+1);
dis=zeros(1,step+1);

for t=1:1
    p=zeros(2*N,step+1);
    v=zeros(2*N,step+1);
    a=zeros(2*N,step+1);
    scal=zeros(1,step+1);

    if N==2
        p(:,1)=[0; 0; 2; 0]*0.625;
        v(:,1)=randi([-100,100],[4,1])/100;
    elseif N==3
        p(:,1)=[0; 0; 2; 0; 4; 0]*0.625;
        v(:,1)=randi([-100,100],[6,1])/100;
    elseif N==5
        p(:, 1)=[0; 0; 2; 0; 4; 0; 1; sqrt(3); 1; -sqrt(3)]*0.625;
        v(:,1)=randi([-100,100],[10,1])/100;
    elseif N==10
        p(:, 1)=[0; 0; 2; 0; 4; 0; 6; 0; 1; sqrt(3); 1; -sqrt(3); 3; sqrt(3); 3; -sqrt(3); 5; sqrt(3); 5; -sqrt(3)]*0.625;
        v(:,1)=randi([-100,100],[20,1])/100;
    elseif N==16
        p(:, 1)=[0; 0; 2; 0; 4; 0; 6; 0; 1; sqrt(3); 1; -sqrt(3); 3; sqrt(3); 3; -sqrt(3); 5; sqrt(3); 5; -sqrt(3); 2; 2*sqrt(3); 2; -2*sqrt(3);
            4; 2*sqrt(3); 4; -2*sqrt(3); 3; 3*sqrt(3); 3; -3*sqrt(3)]*0.625;
        v(:,1)=randi([-100,100],[32,1])/100;
    end
    
    cnt=0;
    for i=1:N
        for j=1:N
            if i~=j
                rp=norm(p(2*j-1:2*j,1)-p(2*i-1:2*i,1));
                if rp<sqrt(d1)
                    dis(1,1)=dis(1,1)+rp;
                    cnt=cnt+1;
                end
                scal(1,1) = scal(1,1) + dot(v(2*i-1:2*i,1),v(2*j-1:2*j,1))/norm(v(2*i-1:2*i,1))/norm(v(2*j-1:2*j,1));
            end
        end
    end
    
    dis(1,1)=dis(1,1)/cnt;

    for i=1:step
        s=0;
        cnt=0;
        d=0;
        for j=1:N
            ta=zeros(2,1);
            f0=zeros(2,1);
            f1=zeros(2,1);
            lambda=0;
            for k=1:N
                if j~=k
                    rp=norm(p(2*j-1:2*j,i)-p(2*k-1:2*k,i));
                    rv=norm(v(2*j-1:2*j,i)-v(2*k-1:2*k,i));
                    if rp<sqrt(d1)
                        d=d+rp;
                        cnt=cnt+1;
                    end
                    if i>j
                        ta=ta + K/(sigma^2+rp^2)^beta * (v(2*k-1:2*k,i)-v(2*j-1:2*j,i));
                        f0=1/(rp^2-d0)^theta*(p(2*j-1:2*j,i)-p(2*k-1:2*k,i));
                        f1=1/(rp^2-d1)^theta*(-p(2*k-1:2*k,i)+p(2*j-1:2*j,i));
                        lambda = lambda+1/N*rv^2;
                    else
                        ta=ta + K/(sigma^2+rp^2)^beta * (v(2*k-1:2*k,i)-v(2*j-1:2*j,i));
                        f0=1/(rp^2-d0)^theta*(p(2*j-1:2*j,i)-p(2*k-1:2*k,i));
                        f1=1/(rp^2-d1)^theta*(-p(2*k-1:2*k,i)+p(2*j-1:2*j,i));
                    end
                end
            end
            a(2*j-1:2*j,i+1)=ta+sqrt(lambda)*(f0-f1);
        end
        v(:,i+1)=v(:,i)+a(:,i+1)*dt;
        p(:,i+1)=p(:,i)+v(:,i)*dt+0.5*a(:,i+1)*dt^2;
        for m=1:N
            for n=1:N
                if m~=n
                    scal(1,i+1) = scal(1,i+1) + dot(v(2*m-1:2*m,i+1),v(2*n-1:2*n,i+1))/norm(v(2*m-1:2*m,i+1))/norm(v(2*n-1:2*n,i+1));
                end
            end
        end
        dis(1,i+1)=d/cnt;
    end
    scal = abs(scal/N/(N-1));
    avg=avg+scal;
end

h=figure;
plot(p(1,:),p(2,:),'o');
hold on;
if N==2
    plot(p(3,:),p(4,:),'o');
elseif N==3
    plot(p(3,:),p(4,:),'o');
    plot(p(5,:),p(6,:),'o');
else
    for i=1:N-1
        plot(p(2*i+1,:),p(2*i+2,:),'o');
    end
end
hold off;

g=figure;
t=1:step+1;
plot(t,dis);
%plot(t,avg/1);

f=figure;
plot(t,n2,'LineWidth',2);
hold on
plot(t,n3,'LineWidth',2)
plot(t,n5,'LineWidth',2)
plot(t,n10,'LineWidth',2)
plot(t,n16,'LineWidth',2)
axis square;
axis([0 step+1 0 1])
xlabel('time steps')
ylabel('\psi_{scal}')
ylabel('average \psi_{scal}')
title('\sigma=1, \beta=0.25')