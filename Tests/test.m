clear all;
close all;
clc;

x = (7:0.01:17);
y1 = trapmf(x,[-inf,-inf,7,11]);

figure();
plot(x,y1);
xlabel('trapmf, P = [1 5 7 11]');
ylim([-0.05 1.05]);

hold on;
y2 = trapmf(x,[7,11,13,17]);
plot(x,y2);

hold on;
y3 = trapmf(x,[13,17,inf,inf]);
plot(x,y3);

hold on;
ytnorm_min1 = min(y1,y2);
ytnorm_min2 = min(y2,y3);
ytnorm_min = ytnorm_min1 + ytnorm_min2;
plot(x,ytnorm_min);

hold on;
[alpha,alpha_i] = max(ytnorm_min);
alpha_mf = y1*0;
alpha_mf(alpha_i) = alpha;
plot(x,alpha_mf);


andmethod = 'no';
A = [y1;y2];
B = AND(A,andmethod);



function out = AND(A,andmethod)

tnorm = ones(1,size(A,2));
if isequal(andmethod, 'min')
    
    for i=1:size(A,1)
        tnorm = min(A(i,:),tnorm);
    end
elseif isequal(andmethod,'product')
    
    for i=1:size(A,1)
        tnorm = A(i,:).*tnorm;
    end
else
    error('incorrect AND method');
end
out = tnorm;
    
end
