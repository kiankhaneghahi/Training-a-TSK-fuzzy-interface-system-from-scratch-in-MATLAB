clear all;
close all;
clc;

x1 = (0:0.01:250);
x2 = (40:0.01:320);
x3 = (0:0.01:100);


trainingmethod= 'Online LSE';
type='tsk1';
numofinputs=3;
numofoutputs=1;
rangeofinputs={x1;x2;x3};
numofmfs=[4;6;2];
membershipfuncstype='gauss';
andmethod='min';
ormethod='max';
defuzzmethod='centroid';
impmethod='min';
aggmethod='sum';

numofdata = 100;

for k = 1 : numofinputs
    
    input_data(k,:) = min(rangeofinputs{k,1}) + (max(rangeofinputs{k,1}) - min(rangeofinputs{k,1}))*rand(1,numofdata,'distributed');
    
end

input_data = gather(input_data);

z = sin(input_data(1,:)) + (sin(input_data(2,:))).^2 + (sin(input_data(3,:))).^3;

io_data = [input_data;z];


K = create_trained_fis(io_data,trainingmethod,type,numofinputs,numofoutputs,rangeofinputs,numofmfs,membershipfuncstype,andmethod,ormethod,defuzzmethod,impmethod,aggmethod);

