%% Mechanism Inversion
% Measurement spheres are iteratively rotated to minimise the fit with
% idealised models
function [output] = focmech(modS,number)

global mechsol
global modFIT
global modnum
modFIT = modS;
modnum = number;

Guess = [0 0 0 0];

options = optimset('MaxFunEvals',1000);

[answer,res,store] = fminsearchbnd(@ellipseMerit,Guess,[0 0 0 -360],[1 1 1 360],options);

global param

output = param;

% Residual calculation
function [res] = ellipseMerit(s)

global modFIT

[store] = ellipseFun(s, modFIT);

test = cell2mat(store(5:8,:))';
test(:,[1,3]) = 1./test(:,[1,3]);
test2(:,1) = test(:,1).*test(:,2);
test2(:,2) = test(:,3).*test(:,4);
test2(:,2) = 1./test2(:,2);

answer = find(test2(:,2)~=0);
answer = answer(test2(answer,2) == min(test2(answer,2)));
res = test2(answer,2);
if isempty(res)
    res = 10000;
end
global param
param = {answer,res,store};

% Model iteration
function [store] = ellipseFun(s, modS)

%store = [];

global modFIT
global models
global modnum
for x = modnum%1:length(models)
    
    fitC = models{x};
    
    if mean(s(:,1:4)) ~= 0
        XYZold = fitC(:,1:3); XYZold = XYZold'; x0=[0 0 0].'; u=[s(:,1) s(:,2)  s(:,3)].'; deg=s(:,4);
        [XYZnew, R, t] = AxelRot(XYZold, deg, u, x0); fitC(:,1:3) = XYZnew';
    end
    
    modS = modFIT;
    IDX = knnsearch(fitC(:,1:3),modS(:,1:3));
    modC = [fitC(IDX,:), modS(:,4)];
    modC(modC(:,4) > 0 & modC(:,5) > 0,6) = 1;
    modC(modC(:,4) < 0 & modC(:,5) < 0,6) = 1;
    modC(1:size(modS,1),7) = 1;
    
    global mechsol
    mechsol = s;
    
    test = []; test2 = [];
    test(1) = sum((modC(:,4) - modC(:,5)).^2);
    test(2) = sum(modC(:,6))/size(modC,1);
    try
        test(3) = sum((modC(1:size(modS,1),4) - (modC(1:size(modS,1),5)./max(abs(modC(1:size(modS,1),5))))).^2);
        test(4) = sum(modC(1:size(modS,1),7))/size(modS,1);
    catch
        test(3) = 10;
        test(4) = 0;
    end
    global modlist
    store{1,x} = modlist{x};
    store{2,x} = [];
    store{3,x} = [s(:,1) s(:,2) s(:,3) s(:,4)];
    store{4,x} = [];
    store{5,x} = test(1);
    store{6,x} = test(2);
    store{7,x} = test(3);
    store{8,x} = test(4);
    store{9,x} = [];
    store{10,x} = [];
    store{11,x} = modC;
    
end