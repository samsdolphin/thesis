clc; clear all
x=0:0.01:10;
y=1.*(x>=0&x<=1)+0;
plot(x,y,'b','linewidth',1.5)
axis([0 1.5 0 1.5])
axis square;
xlabel('||x_i-x_j||^2')
ylabel('a_{ij}(||x_i-x_j||^2)')

x=0:0.01:50;
y=1./(1+x.^2).^0.5;
plot(x,y,'b','linewidth',1.5)
axis([0 50 0 1])
axis square;
xlabel('||x_i-x_j||^2')
ylabel('a_{ij}(||x_i-x_j||^2)')

x=0:0.01:15;
y=1./(x.^2-6).^2;
plot(x,y,'b','linewidth',1.5)
axis([0 10 0 5])
axis square;
xlabel('x')
ylabel('f_0(x)')
line([sqrt(6),sqrt(6)],[0,5],'linestyle','--');