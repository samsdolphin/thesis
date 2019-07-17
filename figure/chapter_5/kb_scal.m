N=2;
beta=1/2;
K=1;
dt=0.5;
sigma=1;
step=1;
p=zeros(2*N,step+1);
v=zeros(2*N,step+1);
a=zeros(2*N,step+1);

if N==2
    p(:,1)=[0; 2; 0; 0];
    v(:,1)=rand(4,1);
elseif N==3
    p(:,1)=[0; 2; -2/sqrt(3); 0; 2/sqrt(3); 0];
    v(:,1)=rand(6,1);
else
    p(:, 1)=[0; 2; -2/sqrt(3); 0; 2/sqrt(3); 0; 0; -2];
    v(:,1)=rand(8,1);
end

for i=1:step
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
end