%clear all;
close all;
N=4;
dt=0.01;
step=300;
d0=1;
d1=2.25;
K=1;
sigma=1;
beta=0.25;
theta=2;

avg_vk=zeros(1,step+1);
avg_cs=zeros(1,step+1);
avg_jiu=zeros(1,step+1);
avg_us=zeros(1,step+1);

dis_vk=zeros(1,step+1);
dis_cs=zeros(1,step+1);
dis_jiu=zeros(1,step+1);
dis_us=zeros(1,step+1);

min_dis_vk=5*ones(1,step);
min_dis_cs=5*ones(1,step);
min_dis_jiu=5*ones(1,step);
min_dis_us=5*ones(1,step);

max_dis_vk=zeros(1,step);
max_dis_cs=zeros(1,step);
max_dis_jiu=zeros(1,step);
max_dis_us=zeros(1,step);

p_vk=zeros(2*N,step+1);
p_cs=zeros(2*N,step+1);
p_jiu=zeros(2*N,step+1);
p_us=zeros(2*N,step+1);

v_vk=zeros(2*N,step+1);
v_cs=zeros(2*N,step+1);
v_jiu=zeros(2*N,step+1);
v_us=zeros(2*N,step+1);

mag_v_vk=zeros(N,step);
mag_v_cs=zeros(N,step);
mag_v_jiu=zeros(N,step);
mag_v_us=zeros(N,step);

a_vk=zeros(2*N,step+1);
a_cs=zeros(2*N,step+1);
a_jiu=zeros(2*N,step+1);
a_us=zeros(2*N,step+1);

scal_vk=zeros(1,step+1);
scal_cs=zeros(1,step+1);
scal_jiu=zeros(1,step+1);
scal_us=zeros(1,step+1);

if N==2
    p_vk(:,1)=[0; 0; 0; 2]*0.625;
    p_cs=p_vk;
    p_jiu=p_vk;
    p_us=p_vk;
    for i=1:N
        ag=rand(1)*pi/2+4/pi;
        v_vk(2*i-1:2*i,1)=sin(ag)*[cos(ag);sin(ag)];
    end
    %v_vk(:,1)=[cos(pi/4);sin(pi/4);cos(3*pi/4);sin(3*pi/4)];
    v_cs=v_vk;
    v_jiu=v_vk;
    v_us=v_vk;
elseif N==3
    p_vk(:,1)=[0; 0; 2; 0; 1; sqrt(3)]*0.625;
    p_cs=p_vk;
    p_jiu=p_vk;
    p_us=p_vk;
    for i=1:N
        ag=rand(1)*pi/2+4/pi;
        v_vk(2*i-1:2*i,1)=[cos(ag);sin(ag)];
    end
    v_vk(:,1)=[cos(1)*cos(pi/4);cos(1)*sin(pi/4);-cos(2)*cos(3*pi/4);-cos(2)*sin(3*pi/4);-cos(3)*cos(-2*pi/4);-cos(3)*sin(-2*pi/4)];
    v_cs=v_vk;
    v_jiu=v_vk;
    v_us=v_vk;
elseif N==4
    p_vk(:,1)=[0; 0; 2; 0; 1; sqrt(3); 1; -sqrt(3)]*0.625;
    p_cs=p_vk;
    p_jiu=p_vk;
    p_us=p_vk;
    for i=1:N
        ag=rand(1)*pi/2+4/pi;
        v_vk(2*i-1:2*i,1)=cos(ag)*[cos(ag);sin(ag)];
    end
    v_vk(:,1)=[cos(1)*cos(pi/4);cos(1)*sin(pi/4);-cos(2)*cos(3*pi/4);-cos(2)*sin(3*pi/4);-cos(3)*cos(3*pi/4);-cos(3)*sin(3*pi/4);-cos(4)*cos(pi/4);-cos(4)*sin(pi/4)];
    v_cs=v_vk;
    v_jiu=v_vk;
    v_us=v_vk;
elseif N==5
    p_us(:, 1)=[0; 0; 2; 0; 4; 0; 1; sqrt(3); 1; -sqrt(3)]*0.625;
    for i=1:N
        ag=rand(1)*pi/2+4/pi;
        v_us(2*i-1:2*i,1)=[cos(ag);sin(ag)];
    end
elseif N==10
    p_vk(:, 1)=[0; 0; 2; 0; 4; 0; 6; 0; 1; sqrt(3); 1; -sqrt(3); 3; sqrt(3); 3; -sqrt(3); 5; sqrt(3); 5; -sqrt(3)]*0.625;
    p_cs=p_vk;
    p_jiu=p_vk;
    p_us=p_vk;
    v_vk(:,1)=randi([-100,100],[20,1])/100;
    v_vk(:,1)=v_start;
    v_cs=v_vk;
    v_jiu=v_vk;
    v_us=v_vk;
elseif N==16
    p_vk(:, 1)=[0; 0; 2; 0; 4; 0; 6; 0; 1; sqrt(3); 1; -sqrt(3); 3; sqrt(3); 3; -sqrt(3); 5; sqrt(3); 5; -sqrt(3); 2; 2*sqrt(3); 2; -2*sqrt(3);
        4; 2*sqrt(3); 4; -2*sqrt(3); 3; 3*sqrt(3); 3; -3*sqrt(3)]*0.625;
    p_cs=p_vk;
    p_jiu=p_vk;
    p_us=p_vk;
    for i=1:N
        ag=rand(1)*2*pi;
        v_vk(2*i-1:2*i,1)=[cos(ag);sin(ag)];
    end
    v_cs=v_vk;
    v_jiu=v_vk;
    v_us=v_vk;
    %v_us(:,1)=randi([-100,100],[32,1])/100;
elseif N==32
    p_vk(:, 1)=[0; 0; 2; 0; 4; 0; 6; 0; 1; sqrt(3); 1; -sqrt(3); 3; sqrt(3); 3; -sqrt(3); 5; sqrt(3); 5; -sqrt(3); 2; 2*sqrt(3); 2; -2*sqrt(3);
        4; 2*sqrt(3); 4; -2*sqrt(3); 3; 3*sqrt(3); 3; -3*sqrt(3); newp];
    p_cs=p_vk;
    p_jiu=p_vk;
    p_us=p_vk;
    for i=1:N
        ag=rand(1)*2*pi;
        v_vk(2*i-1:2*i,1)=[cos(ag);sin(ag)];
    end
    v_cs=v_vk;
    v_jiu=v_vk;
    v_us=v_vk;
    %v_us(:,1)=randi([-100,100],[32,1])/100;
end
cnt_vk=0;
cnt_cs=0;
cnt_jiu=0;
cnt_us=0;
for i=1:N
    for j=1:N
        if i~=j
            rp_vk=norm(p_vk(2*j-1:2*j,1)-p_vk(2*i-1:2*i,1));
            rp_cs=norm(p_cs(2*j-1:2*j,1)-p_cs(2*i-1:2*i,1));
            rp_jiu=norm(p_jiu(2*j-1:2*j,1)-p_jiu(2*i-1:2*i,1));
            rp_us=norm(p_us(2*j-1:2*j,1)-p_us(2*i-1:2*i,1));
            
            if rp_vk<sqrt(4*d1)
                dis_vk(1,1)=dis_vk(1,1)+rp_vk;
                cnt_vk=cnt_vk+1;
            end
            if rp_cs<sqrt(4*d1)
                dis_cs(1,1)=dis_cs(1,1)+rp_cs;
                cnt_cs=cnt_cs+1;
            end
            if rp_jiu<sqrt(4*d1)
                dis_jiu(1,1)=dis_jiu(1,1)+rp_jiu;
                cnt_jiu=cnt_jiu+1;
            end
            if rp_us<sqrt(4*d1)
                dis_us(1,1)=dis_us(1,1)+rp_us;
                cnt_us=cnt_us+1;
            end
            
            scal_us(1,1)=scal_us(1,1) + dot(v_us(2*i-1:2*i,1),v_us(2*j-1:2*j,1))/norm(v_us(2*i-1:2*i,1))/norm(v_us(2*j-1:2*j,1));
            scal_vk(1,1)=scal_us(1,1);
            scal_cs(1,1)=scal_us(1,1);
            scal_jiu(1,1)=scal_us(1,1);
        end
    end
end
dis_vk(1,1)=dis_vk(1,1)/cnt_vk;
dis_cs(1,1)=dis_cs(1,1)/cnt_cs;
dis_jiu(1,1)=dis_jiu(1,1)/cnt_jiu;
dis_us(1,1)=dis_us(1,1)/cnt_us;

for i=1:step
    cnt_vk=0;
    cnt_cs=0;
    cnt_jiu=0;
    cnt_us=0;
    d_vk=0;
    d_cs=0;
    d_jiu=0;
    d_us=0;
    for j=1:N
        ta_vk=zeros(2,1);
        ta_cs=zeros(2,1);
        ta_jiu=zeros(2,1);
        ta_us=zeros(2,1);
        f0_jiu=zeros(2,1);
        f0_us=zeros(2,1);
        f1=zeros(2,1);
        lambda_jiu=0;
        lambda_us=0;
        for k=1:N
            if j~=k
                rp_vk=norm(p_vk(2*j-1:2*j,i)-p_vk(2*k-1:2*k,i));
                rv_vk=norm(v_vk(2*j-1:2*j,i)-v_vk(2*k-1:2*k,i));
                rp_cs=norm(p_cs(2*j-1:2*j,i)-p_cs(2*k-1:2*k,i));
                rv_cs=norm(v_cs(2*j-1:2*j,i)-v_cs(2*k-1:2*k,i));
                rp_jiu=norm(p_jiu(2*j-1:2*j,i)-p_jiu(2*k-1:2*k,i));
                rv_jiu=norm(v_jiu(2*j-1:2*j,i)-v_jiu(2*k-1:2*k,i));
                rp_us=norm(p_us(2*j-1:2*j,i)-p_us(2*k-1:2*k,i));
                rv_us=norm(v_us(2*j-1:2*j,i)-v_us(2*k-1:2*k,i));
                
                if rp_vk<min_dis_vk(1,i)
                    min_dis_vk(1,i)=rp_vk;
                end
                if rp_vk>max_dis_vk(1,i)
                    max_dis_vk(1,i)=rp_vk;
                end
                if rp_cs<min_dis_cs(1,i)
                    min_dis_cs(1,i)=rp_cs;
                end
                if rp_cs>max_dis_cs(1,i)
                    max_dis_cs(1,i)=rp_cs;
                end
                if rp_jiu<min_dis_jiu(1,i)
                    min_dis_jiu(1,i)=rp_jiu;
                end
                if rp_jiu>max_dis_jiu(1,i)
                    max_dis_jiu(1,i)=rp_jiu;
                end
                if rp_us<min_dis_us(1,i)
                    min_dis_us(1,i)=rp_us;
                end
                if rp_us>max_dis_us(1,i)
                    max_dis_us(1,i)=rp_us;
                end
                
                if rp_vk<sqrt(4*d1)
                    d_vk=d_vk+rp_vk;
                    cnt_vk=cnt_vk+1;
                end
                if rp_cs<sqrt(4*d1)
                    d_cs=d_cs+rp_cs;
                    cnt_cs=cnt_cs+1;
                end
                if rp_jiu<sqrt(4*d1)
                    d_jiu=d_jiu+rp_jiu;
                    cnt_jiu=cnt_jiu+1;
                end
                if rp_us<sqrt(4*d1)
                    d_us=d_us+rp_us;
                    cnt_us=cnt_us+1;
                end

                if i>j
                    ta_vk=ta_vk + (v_vk(2*k-1:2*k,i)-v_vk(2*j-1:2*j,i));
                    ta_cs=ta_cs + K/(sigma^2+rp_cs^2)^beta * (v_cs(2*k-1:2*k,i)-v_cs(2*j-1:2*j,i));
                    ta_jiu=ta_jiu + K/(sigma^2+rp_jiu^2)^beta * (v_jiu(2*k-1:2*k,i)-v_jiu(2*j-1:2*j,i));
                    ta_us=ta_us + K/(sigma^2+rp_us^2)^beta * (v_us(2*k-1:2*k,i)-v_us(2*j-1:2*j,i));
                    f0_jiu=1/(rp_jiu^2-d0)^theta*(p_jiu(2*j-1:2*j,i)-p_jiu(2*k-1:2*k,i));
                    f0_us=1/(rp_us^2-d0)^theta*(p_us(2*j-1:2*j,i)-p_us(2*k-1:2*k,i));
                    f1=1/(rp_us^2-d1)^theta*(-p_us(2*k-1:2*k,i)+p_us(2*j-1:2*j,i));
                    lambda_jiu = lambda_jiu+1/N*rv_jiu^2;
                    lambda_us = lambda_us+1/N*rv_us^2;
                else
                    ta_vk=ta_vk + (v_vk(2*k-1:2*k,i)-v_vk(2*j-1:2*j,i));
                    ta_cs=ta_cs + K/(sigma^2+rp_cs^2)^beta * (v_cs(2*k-1:2*k,i)-v_cs(2*j-1:2*j,i));
                    ta_jiu=ta_jiu + K/(sigma^2+rp_jiu^2)^beta * (v_jiu(2*k-1:2*k,i)-v_jiu(2*j-1:2*j,i));
                    ta_us=ta_us + K/(sigma^2+rp_us^2)^beta * (v_us(2*k-1:2*k,i)-v_us(2*j-1:2*j,i));
                    f0_jiu=1/(rp_jiu^2-d0)^theta*(p_jiu(2*j-1:2*j,i)-p_jiu(2*k-1:2*k,i));
                    f0_us=1/(rp_us^2-d0)^theta*(p_us(2*j-1:2*j,i)-p_us(2*k-1:2*k,i));
                    f1=1/(rp_us^2-d1)^theta*(-p_us(2*k-1:2*k,i)+p_us(2*j-1:2*j,i));
                end
            end
        end
        a_vk(2*j-1:2*j,i+1)=ta_vk;
        a_cs(2*j-1:2*j,i+1)=ta_cs;
        a_jiu(2*j-1:2*j,i+1)=ta_jiu+sqrt(lambda_jiu)*f0_jiu;
        a_us(2*j-1:2*j,i+1)=ta_us+sqrt(lambda_us)*(f0_us-f1);
    end
    v_vk(:,i+1)=v_vk(:,i)+a_vk(:,i+1)*dt;
    p_vk(:,i+1)=p_vk(:,i)+v_vk(:,i)*dt+0.5*a_us(:,i+1)*dt^2;
    v_cs(:,i+1)=v_cs(:,i)+a_cs(:,i+1)*dt;
    p_cs(:,i+1)=p_cs(:,i)+v_cs(:,i)*dt+0.5*a_cs(:,i+1)*dt^2;
    v_jiu(:,i+1)=v_jiu(:,i)+a_jiu(:,i+1)*dt;
    p_jiu(:,i+1)=p_jiu(:,i)+v_jiu(:,i)*dt+0.5*a_jiu(:,i+1)*dt^2;
    v_us(:,i+1)=v_us(:,i)+a_us(:,i+1)*dt;
    p_us(:,i+1)=p_us(:,i)+v_us(:,i)*dt+0.5*a_us(:,i+1)*dt^2;
    for m=1:N
        for n=1:N
            if m~=n
                scal_vk(1,i+1) = scal_vk(1,i+1) + dot(v_vk(2*m-1:2*m,i+1),v_vk(2*n-1:2*n,i+1))/norm(v_vk(2*m-1:2*m,i+1))/norm(v_vk(2*n-1:2*n,i+1));
                scal_cs(1,i+1) = scal_cs(1,i+1) + dot(v_cs(2*m-1:2*m,i+1),v_cs(2*n-1:2*n,i+1))/norm(v_cs(2*m-1:2*m,i+1))/norm(v_cs(2*n-1:2*n,i+1));
                scal_jiu(1,i+1) = scal_jiu(1,i+1) + dot(v_jiu(2*m-1:2*m,i+1),v_jiu(2*n-1:2*n,i+1))/norm(v_jiu(2*m-1:2*m,i+1))/norm(v_jiu(2*n-1:2*n,i+1));
                scal_us(1,i+1) = scal_us(1,i+1) + dot(v_us(2*m-1:2*m,i+1),v_us(2*n-1:2*n,i+1))/norm(v_us(2*m-1:2*m,i+1))/norm(v_us(2*n-1:2*n,i+1));
            end
        end
    end
    dis_vk(1,i+1)=d_vk/cnt_vk;
    dis_cs(1,i+1)=d_cs/cnt_cs;
    dis_jiu(1,i+1)=d_jiu/cnt_jiu;
    dis_us(1,i+1)=d_us/cnt_us;
    
    for num=1:N
        mag_v_vk(num,i)=sqrt(v_vk(2*num-1,i)^2+v_vk(2*num,i)^2);
        mag_v_cs(num,i)=sqrt(v_vk(2*num-1,i)^2+v_vk(2*num,i)^2);
        mag_v_jiu(num,i)=sqrt(v_vk(2*num-1,i)^2+v_vk(2*num,i)^2);
        mag_v_us(num,i)=sqrt(v_vk(2*num-1,i)^2+v_vk(2*num,i)^2);
    end
end
scal_vk=abs(scal_vk/N/(N-1));
avg_vk=avg_vk+scal_vk;
scal_cs=abs(scal_cs/N/(N-1));
avg_cs=avg_cs+scal_cs;
scal_jiu=abs(scal_jiu/N/(N-1));
avg_jiu=avg_jiu+scal_jiu;
scal_us=abs(scal_us/N/(N-1));
avg_us=avg_us+scal_us;

g=figure;
plot(p_us(1,:),p_us(2,:),'o');
hold on;
for i=1:N-1
    plot(p_us(2*i+1,:),p_us(2*i+2,:),'o');
end
hold off;
axis equal;
xlabel('x position (m)')
ylabel('y position (m)')
title('N=4')
arrow(p_us(1:2,1),p_us(1:2,1)+v_us(1:2,1),'EdgeColor','b','FaceColor','b')
arrow(p_us(3:4,1),p_us(3:4,1)+v_us(3:4,1),'EdgeColor','r','FaceColor','r')
arrow(p_us(5:6,1),p_us(5:6,1)+v_us(5:6,1),'EdgeColor','y','FaceColor','y')
arrow(p_us(7:8,1),p_us(7:8,1)+v_us(7:8,1),'EdgeColor','y','FaceColor','y')

h=figure;
t=1:step+1;
plot(t,dis_vk,'LineWidth',2);
hold on;
plot(t,dis_cs,'LineWidth',2);
plot(t,dis_jiu,'LineWidth',2);
plot(t,dis_us,'LineWidth',2);
hold off;
axis square;
axis([0 step+1 0 3.5])
xlabel('time steps')
ylabel('average relative distance (m)')
title('N=4')
line([0,step+1],[d0,d0],'linestyle','--');
line([0,step+1],[sqrt(d1),sqrt(d1)],'linestyle','--','color','red');

f=figure;
plot(t,avg_vk,'LineWidth',2);
hold on;
plot(t,avg_cs,'LineWidth',2);
plot(t,avg_jiu,'LineWidth',2);
plot(t,avg_us,'LineWidth',2);
hold off;
axis square;
axis([0 step+1 0 1])
xlabel('time steps')
ylabel('\psi_{angle}')
title('N=4')

mag=figure;
t=1:step;
plot(t,mag_v_vk(1,:),'LineWidth',2);
hold on;
for i=2:N
    plot(t,mag_v_vk(i,:),'LineWidth',2);
end
axis square;
axis([0 step 0 1])
xlabel('time steps')
ylabel('magnitude of velocities')
title('N=2')