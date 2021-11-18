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

% hold on;
% y1y2 = AND([y1;y2],'product');
% plot(x,y1y2);
% 
% hold on;
% y2y3 = AND([y2;y3],'product');
% plot(x,y2y3);

hold on;
y1y2 = OR([y1;y2],'max');
plot(x,y1y2);

hold on;
y2y3 = OR([y2;y3],'max');
plot(x,y2y3)


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

function out = OR(A,ormethod)

snorm = zeros(1,size(A,2));
if isequal(ormethod, 'max')
    
    for i=1:size(A,1)
        snorm = max(A(i,:),snorm);
    end
elseif isequal(ormethod,'sum')
    
    for i=1:size(A,1)
        snorm = A(i,:)+snorm;
    end
else
    error('incorrect OR method');
end
out = snorm;
    
end