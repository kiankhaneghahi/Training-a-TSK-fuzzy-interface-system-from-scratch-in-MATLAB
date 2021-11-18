clear all;
close all;
clc;


numofinputs = 3;
numofmfs = [6;4;2];

rulebase = Rule_Base(numofinputs,numofmfs);

for rulenum = 1 : size(rulebase,2)
    
    input_disp = [];
    
    for k = 1 : size(rulebase,1)
        
        input_disp = [input_disp,' x',num2str(k),' is MF',num2str(k),num2str(rulebase(k,rulenum))];
        if k ~= size(rulebase,1)
            
            input_disp = [input_disp,' and '];
           
        end
        
    end
    
    disp(['If ',input_disp,' Then ']);
    
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