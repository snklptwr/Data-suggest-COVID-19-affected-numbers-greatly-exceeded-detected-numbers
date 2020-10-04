function z = ContinuumMainFunction(X)
global nDays data_cal data country ys sigma cf gamma1 tolN
% to reproduce Table 3 in the paper, enter
% X=[-2.00188754602483,1.62595165533504,1.95047276782753,2.29892201682886,-1.61792666323166] for Italy
% X=[-1.89522721309271,1.67228015713143,1.97912619743166,-0.153459534071443,-2.54704622412673] for Germany
% X=[-2.13319452757009,1.58964187921160,1.71034216222796,-1.46787516089771,-0.740241543278229] for UK
% X=[-2.94303913317586,1.78072240689756,1.31624278957644,1.00759547374750,-4.00423791046035] for Spain
sigma=3; gamma1=0.07; tolN=8; country='Italy';
[data,ys,cf]=countryData(country); 
data=data/ys/1e8; nDays=numel(data); 
data_cal=ddeOptimContinuum(X);
e=data_cal-data(cf:end);
z=norm(e)/norm(data(cf:end));
z=100*z;
end