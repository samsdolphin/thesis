%clear all;
close all;
N=3;
dt=0.01;
step=300;
d0=0.5;
d1=9;
K=1;
sigma=1;
beta=0.25;
theta=2;

avg=zeros(1,step+1);
dis=zeros(1,step+1);
p=zeros(2*N,step+1);
v=zeros(2*N,step+1);
a=zeros(2*N,step+1);
scal=zeros(1,step+1);

if N==2
    p(:,1)=[0; 0; 0; 2]*0.625;
    %v=[cos(pi/4);sin(pi/4);cos(3*pi/4);sin(3*pi/4)];
    for i=1:N
        ag=rand(1)*pi/2+4/pi;
        v(2*i-1:2*i,1)=[cos(ag);sin(ag)];
    end
elseif N==3
    p(:,1)=[0; 0; 2; 0; 1; sqrt(3)]*0.625;
    %p(:,1)=[0; 0; 2; 0; 4; 0]*0.625;    
    for i=1:N
        ag=rand(1)*pi/2+4/pi;
        v(2*i-1:2*i,1)=[cos(ag);sin(ag)];
    end
    %v(:,1)=[cos(pi/4);sin(pi/4);3*cos(2*pi/4);3*sin(2*pi/4);cos(3*pi/4);sin(3*pi/4)];
elseif N==4
    p(:,1)=[0; 0; 2; 0; 1; sqrt(3); 1; -sqrt(3)]*0.625;
    %v(:,1)=[cos(pi/4);sin(pi/4);cos(3*pi/4);sin(3*pi/4);cos(2*pi/4);sin(2*pi/4);cos(pi/4);sin(pi/4)];
    for i=1:N
        ag=rand(1)*pi/2+4/pi;
        v(2*i-1:2*i,1)=[cos(ag);sin(ag)];
    end
elseif N==5
    p(:, 1)=[0; 0; 2; 0; 4; 0; 1; sqrt(3); 1; -sqrt(3)]*0.625;
    for i=1:N
        ag=rand(1)*pi/2+4/pi;
        v(2*i-1:2*i,1)=[cos(ag);sin(ag)];
    end
elseif N==10
    p(:, 1)=[0; 0; 2; 0; 4; 0; 6; 0; 1; sqrt(3); 1; -sqrt(3); 3; sqrt(3); 3; -sqrt(3); 5; sqrt(3); 5; -sqrt(3)]*0.625;
    for i=1:N
        ag=rand(1)*pi/2+4/pi;
        v(2*i-1:2*i,1)=[cos(ag);sin(ag)];
    end
elseif N==16
    p(:, 1)=[0; 0; 2; 0; 4; 0; 6; 0; 1; sqrt(3); 1; -sqrt(3); 3; sqrt(3); 3; -sqrt(3); 5; sqrt(3); 5; -sqrt(3); 2; 2*sqrt(3); 2; -2*sqrt(3);
        4; 2*sqrt(3); 4; -2*sqrt(3); 3; 3*sqrt(3); 3; -3*sqrt(3)]*0.625;
    for i=1:N
        ag=rand(1)*pi/2+4/pi;
        v(2*i-1:2*i,1)=[cos(ag);sin(ag)];
    end
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
            scal(1,1)=scal(1,1) + dot(v(2*i-1:2*i,1),v(2*j-1:2*j,1))/norm(v(2*i-1:2*i,1))/norm(v(2*j-1:2*j,1));
        end
    end
end

dis(1,1)=dis(1,1)/cnt;

for i=1:step
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

scal=abs(scal/N/(N-1));
avg=avg+scal;

g=figure;
plot(p(1,:),p(2,:),'o');
hold on;
for i=1:N-1
    plot(p(2*i+1,:),p(2*i+2,:),'o');
end
hold off;
axis equal;
arrow(p(1:2,1),p(1:2,1)+v(1:2,1),'EdgeColor','b','FaceColor','b')
arrow(p(3:4,1),p(3:4,1)+v(3:4,1),'EdgeColor','r','FaceColor','r')
arrow(p(5:6,1),p(5:6,1)+v(5:6,1),'EdgeColor','y','FaceColor','y')
line([p(1,end),p(3,end)],[p(2,end),p(4,end)],'linestyle','--');
line([p(1,end),p(5,end)],[p(2,end),p(6,end)],'linestyle','--');
line([p(3,end),p(5,end)],[p(4,end),p(6,end)],'linestyle','--');

h=figure;
t=1:step+1;
plot(t,dis,'LineWidth',2);
axis square;

f=figure;
plot(t,avg,'LineWidth',2);

axis square;
axis([0 step+1 0 1])
xlabel('time steps')
ylabel('average \psi_{scal}')
title('N=2')