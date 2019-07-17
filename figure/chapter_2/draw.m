% Vicsek's model aij
x=0:0.01:1.5;
y=1*(x>=0&x<=1);
plot(x,y,'b','linewidth',1.5)
axis equal
axis([0 1.5 0 1.5])

% Cucker Smale's model aij
x=0:0.1:50;
fplot(@(x) 1/(1+x^2)^0.45,'b','linewidth',1.5)
axis([0 50 0 1])
axis square

% Dong Cucker's model f(r)
x=1.1:0.1:10;
fplot(@(x) 1/(x-1)^2,'b','linewidth',1.5)
axis([0 15 0 5])
axis square

f = @(x,y) (x-0).^2 +(y-0).^2;
g = @(x,y) 1./((x-3).^2+(y-1.5).^2+0.05)+1./((x-1.5).^2+(y-3).^2+0.05)+1./((x-1.5).^2+(y-1.5).^2+0.05)+1./((x-3).^2+(y-3).^2+0.05);
[X,Y] = meshgrid(0:0.05:5);
figure(1)
meshc(X, Y, f(X,Y)+g(X,Y))
grid on