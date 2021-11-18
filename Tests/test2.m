clear all;
close all;
clc;

x = (0:0.01:250);
x1 = (0:0.01:100);
numofinputs = 3;
rangeofinputs = {x;x;x1};
numofmfs = [8;5;3];
membershipfuncstype = 'bell';


inputs = create_inputs(numofinputs,rangeofinputs,numofmfs,membershipfuncstype);

figure(1);
plot(x,inputs{1,1});
xlabel('input1');
ylim([-0.05 1.05]);


figure(2);
plot(x,inputs{2,1});
xlabel('input2');
ylim([-0.05 1.05]);

figure(3);
plot(x1,inputs{3,1});
xlabel('input3');
ylim([-0.05 1.05]);







function inputsout = create_inputs(numofinputs,rangeofinputs,numofmfs,membershipfuncstype)
inputs = {};

for k = 1:numofinputs
    x = rangeofinputs{k,:};
    xmin = min(x);
    xmax = max(x);
    R = xmax -xmin;
   
    
        
    if isequal(membershipfuncstype,'trapezoid')

        delta = R/(2*(numofmfs(k,:)-1));
        mf = trapmf(x,[-inf,-inf,xmin+delta/2,xmin+3*delta/2]);
        a1 = xmin + delta/2; a2 = xmin + 3*delta/2; a3 = xmin + 5*delta/2; a4 = xmin + 7*delta/2;
        for i = 2:numofmfs(k,:)-1
            
            mf_i = trapmf(x,[a1,a2,a3,a4]);
            mf = [mf;mf_i];
            
            a1 = a1 + 2*delta;
            a2 = a2 + 2*delta;
            a3 = a3 + 2*delta;
            a4 = a4 + 2*delta;
            
        end
        mf_final = trapmf(x,[xmax-3*delta/2,xmax-delta/2,inf,inf]);
        mf = [mf;mf_final];

    elseif isequal(membershipfuncstype,'triangle')
        
        delta = R/(numofmfs(k,:)-1);
        mf = trimf(x,[-inf,xmin,xmin+delta]);
        a1 = xmin; a2 = xmin + delta; a3 = xmin + 2*delta;
        for i = 2:numofmfs(k,:)-1
            
            mf_i = trimf(x,[a1,a2,a3]);
            mf = [mf;mf_i];
            
            a1 = a1 + delta;
            a2 = a2 + delta;
            a3 = a3 + delta;
            
        end
        mf_final = trimf(x,[xmax-delta,xmax,inf]);
        mf = [mf;mf_final];

    elseif isequal(membershipfuncstype,'bell')
        
        delta = R/(2*(numofmfs(k,:)-1));
        a = delta*0.866;
        b = 5.45367;
        
        mf = gbellmf(x,[a,b,xmin]);
        c = xmin + 2*delta;
        for i = 2:numofmfs(k,:)-1
            
            mf_i = gbellmf(x,[a,b,c]);
            mf = [mf;mf_i];
            
            c = c + 2*delta;
            
        end
        mf_final = gbellmf(x,[a,b,xmax]);
        mf = [mf;mf_final];

    elseif isequal(membershipfuncstype,'gauss')
        
        delta = R/(numofmfs(k,:)-1);
        sigma = delta/3.5;
        mf = gaussmf(x,[sigma,xmin]);
        c = xmin + delta;
        for i = 2:numofmfs(k,:)-1
            
            mf_i = gaussmf(x,[sigma,c]);
            mf = [mf;mf_i];
            
            c = c + delta;
            
        end
        mf_final = gaussmf(x,[sigma,xmax]);
        mf = [mf;mf_final];

    else
         error('Incorrect membership function type');
    end
    
    inputs(k,1)= {mf};
end

inputsout = inputs;

end