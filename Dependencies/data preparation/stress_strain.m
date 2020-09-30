%% Load Mechanical Stress-Strain Data
% Stress-strain data are loaded, corrected for machine stiffness and
% time-shifted so that the experiment begins at zero time.

%% Load data
% Strain is stored in deform.txt as two column matrix of time and strain
% Differential stress is stored in stress.txt

load deform.txt
load stress.txt
try
    load force.txt
end




%% Machine stiffness correction
% Stiffness correction is made according to the PhD thesis of Marco Fazio
try
    deform(:,2) = deform(:,2) - (stress(:,2)./130000);
catch
    stress(end,:) = [];
    deform(:,2) = deform(:,2) - (stress(:,2)./130000);
end
%% Data windowing
% Data are windowed between the first instance of stress exceeding 1 MPa
% and strain reaching its maximum value. Data can be windowed in a number
% of ways shown by ind2

ind = find(stress(:,2) >= 5);
ind = ind(1);
if isempty(ind) == 1
    ind = 1;
end

% End at UCS

%ind2 = find(stress(:,2) >= max(stress(:,2))); 
% ind2 = find(deform(:,2) >= round((deform(ind2(1),2))/0.01)*0.01);
% ind2 = find(deform(:,2) >= deform(ind2(1),2)+0.1);

% End just after failure

%ind2 = find(diff(deform(:,2)) == max(diff(deform(:,2))))+2; 
%ind2 = find(deform(:,2) >= (round((deform(ind2(1),2))/0.1))*0.1);

% End at maximum strain

ind2 = find(deform(ind(1):end,2) >= max(deform(ind(1):end,2))-0.1)+ind(1);

deform = deform(ind(1):ind2(1),:);
stress = stress(ind(1):ind2(1),:);
try
    force = force(ind(1):ind2(1),:);
end

deform(:,2) = deform(:,2) - min(deform(:,2));
stress(:,2) = stress(:,2) - min(stress(:,2));

%% Time correction

start = min(deform(1,1));
deform(:,1) = deform(:,1)- start;
stress(:,1) = deform(:,1);

try
    force(:,1) = deform(:,1);
end