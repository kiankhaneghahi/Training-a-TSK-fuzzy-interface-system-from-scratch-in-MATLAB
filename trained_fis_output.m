function out = trained_fis_output()



end




function P_i_out = P_i(rangeofinputs,inputs,input_data,rule_base,numofinputs,n,m,M,aggmethod,andmethod)

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