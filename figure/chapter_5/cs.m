close all;
N=2;
dt=0.5;
step=50;
K=1;
sigma=1;
beta=0.25;
loop=1;
avg=zeros(1,step+1);
for t=1:loop
    p=zeros(2*N,step+1);
    v=zeros(2*N,step+1);
    a=zeros(2*N,step+1);
    scal=zeros(1,step+1);

    if N==2
        p(:,1)=[0; 2; 0; 0];
        v(:,1)=randi([-100,100],[4,1])/100;
    elseif N==3
        p(:,1)=[0; 2; -2/sqrt(3); 0; 2/sqrt(3); 0];
        v(:,1)=randi([-100,100],[6,1])/100;
    else
        p(:, 1)=[0; 2; -2/sqrt(3); 0; 2/sqrt(3); 0; 0; -2];
        v(:,1)=randi([-100,100],[8,1])/100;
    end

    for i=1:N
        for j=1:N
            if i~=j
                scal(1,1) = scal(1,1) + dot(v(2*i-1:2*i,1),v(2*j-1:2*j,1))/norm(v(2*i-1:2*i,1))/norm(v(2*j-1:2*j,1));
            end
        end
    end

    for i=1:step
        s=0;
        for j=1:N
            ta=zeros(2,1);
            for k=1:N
                if j~=k
                    rp=norm(p(2*j-1:2*j,i)-p(2*k-1:2*k,i));
                    ta=ta + K/(sigma^2+rp^2)^beta * (v(2*k-1:2*k,i)-v(2*j-1:2*j,i));
                end
            end
            a(2*j-1:2*j,i+1)=ta;
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
    plot(p(3,:),p(4,:),'o');
    plot(p(5,:),p(6,:),'o');
    plot(p(7,:),p(8,:),'o');
end
hold off;

g=figure;
t=1:step+1;
plot(t,avg/loop);
axis square;
axis([0 step+1 0 1])
xlabel('time steps')
ylabel('\psi_{scal}')
ylabel('average \psi_{scal}')
title('\sigma=0.5')