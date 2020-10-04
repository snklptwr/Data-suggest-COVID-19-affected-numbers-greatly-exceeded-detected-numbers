function [data, ys, cf] = countryData(country)

if strcmp(country,'Italy')
    load(strcat(country,'_TotalCases.mat'),strcat(country,'_TotalCases'));
    data=Italy_TotalCases; nDays=numel(data); nDArray=1:nDays;
    ys=0.6; cf=min(nDArray(data>data(end)*0.01)); 
elseif strcmp(country,'Germany')
    load(strcat(country,'_TotalCases.mat'),strcat(country,'_TotalCases'));
    data=Germany_TotalCases; nDays=numel(data); nDArray=1:nDays;
    ys=0.837; cf=min(nDArray(data>data(end)*0.01));  
elseif strcmp(country,'Spain')
    load(strcat(country,'_TotalCases.mat'),strcat(country,'_TotalCases'));
    data=Spain_TotalCases; nDays=numel(data); nDArray=1:nDays;
    ys=0.469; cf=min(nDArray(data>data(end)*0.01)); 
elseif strcmp(country,'UK')
    load(strcat(country,'_TotalCases.mat'),strcat(country,'_TotalCases'));
    data=UK_TotalCases; nDays=numel(data); nDArray=1:nDays; 
    ys=0.666; cf=min(nDArray(data>data(end)*0.01));
end
    