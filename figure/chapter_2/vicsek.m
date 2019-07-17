close all;
n = 5;
x = 2*rand(2,n);
v = abs(rand(2,n));
%v=0.5*ones(2,n)+0.1*abs(rand(2,n));
a = zeros(2,n);
T = 5;
dt = 0.1;
X=zeros(2*n,T/dt);

h=figure;
for k=1:T/dt
    a = zeros(2,n);
    for i=1:n
        for j=1:n
            if i~=j && norm(x(:,i)-x(:,j))<=1
                a(:,i) = a(:,i) + v(:,j)-v(:,i);
            end
        end
    end
    v = v + a*dt;
    x = x + v*dt + 0.5*a*dt.^2;
    for i=1:n
        X(i,k) = x(1,i);
        X(i+n,k) = x(2,i);
    end
end
hold on;
plot(X(1,:),X(1+n,:),'c^');
plot(X(2,:),X(2+n,:),'m^');
plot(X(3,:),X(3+n,:),'g^');
plot(X(4,:),X(4+n,:),'b^');
plot(X(5,:),X(5+n,:),'r^');
axis equal;
axis([0 5 0 5])
rectangle('Position',[0 0 1.5 1.5])