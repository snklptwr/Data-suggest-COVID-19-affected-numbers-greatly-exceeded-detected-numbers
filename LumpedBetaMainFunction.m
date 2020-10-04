function z = LumpedBetaMainFunction(X)
global nDays data_cal data country ys sigma cf gamma1 tolN
% to reproduce Table 3 in the paper, enter
% X=[-1.70077176705963,2.11973581888038,0.395719351626882,-0.721260697939164] for Italy
% X=[-1.56197393248550,2.23182335036030,-0.119353967471291,-0.651883083842769] for Germany
% X=[-1.81013042141519,2.07079833806585,0.385162980548601,-0.739544469042364] for UK
% X= [-1.72308172558173,1.98492630495309,-0.343679650958787,-0.302993722634888] for Spain
sigma=3; gamma1=0.07; tolN=8; country='Italy';
[data,ys,cf]=countryData(country); 
data=data/ys/1e8; nDays=numel(data); 
data_cal=ddeOptimLumped(X);
e=data_cal-data(cf:end);
z=norm(e)/norm(data(cf:end));
z=100*z;
end