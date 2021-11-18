

function [TrainingParameters,RuleBase,inputs] = create_trained_fis(io_data,trainingmethod,type,numofinputs,numofoutputs,rangeofinputs,numofmfs,membershipfuncstype,andmethod,ormethod,defuzzmethod,impmethod,aggmethod)


% create_trained_fis(io_data,trainigmethod,type,numofinputs,numofoutputs,rangeofinputs,numofmfs,membershipfuncstype,andmethod,ormethod,defuzzmethod,impmethod,aggmethod)
%
%
% trainingmethod = 'Offline LSE' or 'Online LSE' or 'GD'
% type = 'mamdani' or 'tsk0' or 'tsk1'
% numofinputs = 1 or 2 or ...
% numofoutputs = 1 or 2 or ...
% rangeofinputs = {range1,range2,...}
% numofmfs = [2 or 3 or ... , 2 or 3 or ..., ...]
% membershipfuncstype = 'guass' or 'bell' or 'triangle' or 'trapezoid'
% andmethod = 'min' or 'product'
% ormethod = 'max' or 'sum'
% defuzzmethod = 'centroid' or 'bisector' or 'mom' or 'lom' or 'sum' - ('mom' — Mean of the values for which the output fuzzy set is maximum | 'lom' — Largest value for which the output fuzzy set is maximum | 'som' — Smallest value for which the output fuzzy set is maximum)
% impmethod = 'min' or 'product'
% aggmethod = 'max' or 'sum'

numofrules = 1;

for i=1:numofinputs
    
    numofrules = numofmfs(i) * numofrules;
    
end

inputs = create_inputs(numofinputs,rangeofinputs,numofmfs,membershipfuncstype);

x = {};
for k = 1:numofinputs
    
    x_i = rangeofinputs{k,:};
    x = {x;x_i};
    figure(k);
    plot(x_i,inputs{k,1});
    xlabel(['input ',num2str(k)]);
    ylabel(['\mu(',num2str(k),')']);
    ylim([-0.05 1.05]);

end

rule_base = Rule_Base(numofinputs,numofmfs);
RuleBase = rule_base;

if isequal(type,'tsk0')
    
    K = trained_tsk0_parameters(io_data,rule_base,inputs,trainingmethod,numofinputs,numofoutputs,numofrules,numofmfs,rangeofinputs,andmethod,aggmethod);
    print_rulebase(rule_base,K,type);
    TrainingParameters = K;

elseif isequal(type,'tsk1')
    
    K = trained_tsk1_parameters(io_data,rule_base,inputs,trainingmethod,numofinputs,numofoutputs,numofrules,numofmfs,rangeofinputs,andmethod,aggmethod);
    print_rulebase(rule_base,K,type);
    TrainingParameters = K;
    
elseif isequal(type,'mamdani')
    
    error('Mamdani FIS type is still under development.');
    
else
    
    error('Incorrect FIS type.');
    
end
    

end




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




function rule_base = Rule_Base(numofinputs,numofmfs)

numofmfcombinations = 1;

for k = 1 : numofinputs
    
    numofmfcombinations = numofmfcombinations * numofmfs(k,:);
    
end

rulebase = zeros(numofinputs,numofmfcombinations);

for k = 1 : numofinputs
    
    
    other_combinations_num = 1;

    for k_other = k+1 : numofinputs

        other_combinations_num = other_combinations_num * numofmfs(k_other);

    end
    
    counter = 1;
    
    while(counter <= numofmfcombinations)

        for mfnumber = 1 : numofmfs(k,:)

            for combination = 1 : other_combinations_num

                rulebase(k,counter) = mfnumber;
                counter = counter + 1;

            end

        end
    
    end
        
    
end

rule_base = rulebase;

end




function Theta = trained_tsk0_parameters(io_data,rule_base,inputs,trainingmethod,numofinputs,numofoutputs,numofrules,numofmfs,rangeofinputs,andmethod,aggmethod)


for k = 1 : numofinputs
    
    input_data(k,:) = io_data(k,:);
    
end

for k = 1 : numofoutputs
    
    output_data(k,:) = io_data(k+numofinputs,:);
    
end

N = size(input_data,2);
M = numofrules;

Phi = zeros(N,M);

for n = 1 : N

    for m = 1 : M
        
        Phi(n,m) = P_i(rangeofinputs,inputs,input_data,rule_base,numofinputs,numofmfs,n,m,M,aggmethod,andmethod);

    end

end

if isequal(trainingmethod,'Offline LSE')
    
    theta = inv(Phi'*Phi)*(Phi'*output_data');
    
elseif isequal(trainingmethod,'Online LSE')
    
    total_parameters_num = numofrules;
    theta = Online_LSE(Phi,total_parameters_num,output_data);
    
elseif isequal(trainingmethod,'GD')
    
    error('Gradient Descent training method is still under development.');
    
else
    
    error('Incorrect training method.');

end

Theta = theta;

end




function Theta = trained_tsk1_parameters(io_data,rule_base,inputs,trainingmethod,numofinputs,numofoutputs,numofrules,numofmfs,rangeofinputs,andmethod,aggmethod)


for k = 1 : numofinputs
    
    input_data(k,:) = io_data(k,:);
    
end

for k = 1 : numofoutputs
    
    output_data(k,:) = io_data(k+numofinputs,:);
    
end

N = size(input_data,2);
M = numofrules;

Phi = zeros(N,M*(numofinputs+1));

for n = 1 : N

    for m = 1 : M

        for k = 1 : numofinputs
            
            Phi(n,(numofinputs+1)*(m-1)+k) = P_i(rangeofinputs,inputs,input_data,rule_base,numofinputs,numofmfs,n,m,M,aggmethod,andmethod)*input_data(k,n);
            
        end
        
        Phi(n,(numofinputs+1)*(m-1)+numofinputs+1) = P_i(rangeofinputs,inputs,input_data,rule_base,numofinputs,numofmfs,n,m,M,aggmethod,andmethod);

    end

end

if isequal(trainingmethod,'Offline LSE')
    
    theta = inv(Phi'*Phi)*(Phi'*output_data');
    
elseif isequal(trainingmethod,'Online LSE')
    
    total_parameters_num = (numofinputs + 1)*numofrules;
    theta = Online_LSE(Phi,total_parameters_num,output_data);
    
elseif isequal(trainingmethod,'GD')
    
    error('Gradient Descent training method is still under development.');
    
else
    
    error('Incorrect training method.');

end

Theta = theta;

end




function P_i_out = P_i(rangeofinputs,inputs,input_data,rule_base,numofinputs,numofmfs,n,m,M,aggmethod,andmethod)

inputs_data_n = input_data(:,n);
rule_m = rule_base(:,m);

denom = 1;

for m = 1 : M
    
    alpha = 1;
    
    for k = 1 : numofinputs
        
        x_k = rangeofinputs{k,1};
        inputs_k = inputs{k,1};
        stepsize = (max(x_k) - min(x_k))/(size(x_k,2)-1);
        impnum = round((inputs_data_n(k,:)-min(x_k)+stepsize)/stepsize);
        alpha = AND([alpha;inputs_k(rule_base(k,m),impnum)],andmethod);
        
    end
    
    denom = AGG([denom; alpha],aggmethod);

end

alpha_i = 1;

for k = 1 : numofinputs
        
        x_k = rangeofinputs{k,1};
        inputs_k = inputs{k,1};
        stepsize = (max(x_k) - min(x_k))/(size(x_k,2)-1);
        impnum = round((inputs_data_n(k,:)-min(x_k)+stepsize)/stepsize);
        alpha_i = AND([alpha_i;inputs_k(rule_m(k,:),impnum)],andmethod);
        
end

numerator = alpha_i;

P_i_out = numerator/denom;

end




function theta = Online_LSE(A,total_parameters_num,output_data)

y = output_data';      % Output Data
m = total_parameters_num; 

n = size(output_data,2);        % Data Size
alpha = 1000000;
P0 = alpha*eye(m);
theta0 = zeros(m,1);

Pk_old = P0;
theta_old = theta0;

for k = 1:n
    
    Pk_new = Pk_old - (Pk_old*A(k,:)'*A(k,:)*Pk_old)/(1+A(k,:)*Pk_old*A(k,:)');
    theta_new = theta_old + Pk_new*A(k,:)'*(y(k)-A(k,:)*theta_old);
    theta_old = theta_new;
    Pk_old = Pk_new;
    
end

theta = theta_new;

end




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




function out = IMP(A,impmethod)

tnorm = ones(1,size(A,2));
if isequal(impmethod, 'min')
    
    for i=1:size(A,1)
        tnorm = min(A(i,:),tnorm);
    end
elseif isequal(impmethod,'product')
    
    for i=1:size(A,1)
        tnorm = A(i,:).*tnorm;
    end
else
    error('incorrect IMP method');
end
out = tnorm;
    
end




function out = AGG(A,aggmethod)

snorm = zeros(1,size(A,2));
if isequal(aggmethod, 'max')
    
    for i=1:size(A,1)
        snorm = max(A(i,:),snorm);
    end
elseif isequal(aggmethod,'sum')
    
    for i=1:size(A,1)
        snorm = A(i,:)+snorm;
    end
else
    error('incorrect AGG method');
end
out = snorm;
    
end




function print_rulebase(rule_base,K,tsk_type)

% tsk_type = 'tsk0' or 'tsk1'

numofrules = size(rule_base,2);
numofinputs = size(rule_base,1);

for rulenum = 1 : numofrules
    
    input_disp = [];
    output_disp = [];
    
    for k = 1 : numofinputs
        
        input_disp = [input_disp,'x',num2str(k),' is MF',num2str(k),num2str(rule_base(k,rulenum))];
        
        if k ~= numofinputs
            
            input_disp = [input_disp,' and '];
           
        end
        
    end
    if isequal(tsk_type,'tsk0')
        
        for z = 1 : size(K,2)

            output_disp = [output_disp,'Z',num2str(z),' = ',];

            output_disp = [output_disp,num2str(K(rulenum,z))];

            if z ~= size(K,2)

                output_disp = [output_disp,' and '];

            end

        end
        
    elseif isequal(tsk_type,'tsk1')
        
        for z = 1 : size(K,2)

            output_disp = [output_disp,'Z',num2str(z),' = ',];

            for k = 1 : numofinputs

                output_disp = [output_disp,num2str(K((rulenum-1)*(numofinputs+1)+k,z)),'*x',num2str(k)];

                output_disp = [output_disp,' + '];

            end

            output_disp = [output_disp,num2str(K((rulenum-1)*(numofinputs+1)+numofinputs+1,z))];

            if z ~= size(K,2)

                output_disp = [output_disp,' and '];

            end

        end
        
    else
        
        error('Incorrect TSK type');
        
    end
    
    disp(['Rule ',num2str(rulenum),' : If ',input_disp,' Then ',output_disp]);
    
end

end



