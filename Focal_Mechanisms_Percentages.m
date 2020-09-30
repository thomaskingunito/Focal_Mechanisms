%% Mechanism Probability Density Plots
% This code plots the focal mechanisms as probability densities against
% time or strain

%% Version
% Version 1.0, 27th January 2019. Thomas King
%   - First Version

%% Parameter customisation
% Below are the suggested parameters to be modified. I don't recommend
% changing any of the code outside of these parameters.

clear all;% close all

% Plotting Parameters
pressure = '20 MPa'; % This is the title of the plot
TorS = 2; % time or strain 1 or 2
averagepolarity = 0; % 1 on 0 off
ampthresh = 0.05; % Seperate mechanisms with an amplitude threshold
saving = 0;

% Plotting colours. Have the same number of colours as mechanisms
C = brighten(parula(3),.25);
C1 = C(1,:);
C2 = brighten(C(2,:),-.1);
C3 = C(3,:);
C = [C1;C2;C3];

% Smoothing Parameters
pdfsmooth = 40; % PDF plot smoothing
nEvents = 10; % Event windowing

% Time Corrections
timecorr = 0;

%% Compile data

% load mechanical data
stress_strain

% Load and order data
load eventdatamech_ml_residual.mat
%%
[~,order] = sort(cell2mat(eventdata(:,2)));
eventdata = eventdata(order,:);

% Mechanism list
modlist = {'fitCLVD','fitDCQ','fitMM'};
%modlist = {'fitMM','fitDCQ','fitCLVD'};

% Compile mechanism data
ind = [];
load focalmechmodel.mat fitMM
for i = 1:size(eventdata,1)
    
    % Skip unsolved mechanisms
    if isempty(eventdata{i,14}) == 1 || isempty(eventdata{i,2}) == 1
        continue
    end
    
    % Time correction
    eventdata{i,2} = eventdata{i,2} - start - timecorr;
    eventtime(i) = eventdata{i,2};
    
    % Load data for event
    store = eventdata{i,13}; % Fitting parameters
    amp = eventdata{i,8}; % Polarity amplitude
    pol = eventdata{i,9}; % Polarity direction
    avepol = mean(pol(pol~=0)); % Average Polarity
    
    % Maximum amplitude of event
    csig = eventdata{i,6}; csig = csig(:,rms(csig)==max(rms(csig)));
    csig = log(max(abs(csig)));
    
    % Fitting
    test = cell2mat(store(5:8,:))';
    test(:,[1,3]) = 1./test(:,[1,3]);
    test2(:,1) = test(:,1).*test(:,2);
    test2(:,2) = test(:,3).*test(:,4);
    [~,order] = sort(test2(:,2),'descend');
    eventdata{i,16} = test(order(1),3);
    eventdata{i,17} = test(order(1),4);
    eventdata{i,18} = test2(order(1),2);
    fitvalue(i) = test2(order(1),2);
    
%     if fitvalue(i) < 2
%         continue
%     end
    
    % Amplitude data
    aT2(i,find(ismember(modlist,eventdata(i,14)) == 1)) = csig;
    
    aT3(i) = csig;
    
    % Average Polarity Fitting
    if averagepolarity == 1
        if avepol < -0.25
            eventdata{i,14} = modlist{4};
        elseif avepol > 0.25
            eventdata{i,14} = modlist{1};
        elseif avepol >= -0.25 && avepol <= 0.25
            eventdata{i,14} = modlist{3};
        end
    end
    
    % Removes skipped data
    ind = [ind,i];
    eventdata{i,15} = length(ind);
    
    % Converts event time to strain value
    straintime(i,1) = mean(deform(abs(deform(:,1) - eventdata{i,2}) == min(abs(deform(:,1) - eventdata{i,2})),2));
    
end

% Cropping
eventdata = eventdata(ind,:);
straintime = straintime(ind);
aT2 = aT2(ind,:);
fitvalue = fitvalue(ind);
aT3 = aT3(ind);
eventtime =  eventtime(ind);

%% Plotting

% Use this to choose specific events
indE = [1:1:size(eventdata,1)];

% Amplitude thresholding
mechsep = 2.*[1:1:size(modlist,2)]-1; ls = []; ls2 = [];
for j = 1:length(indE) 
    % Compile event data
    ls(j,1) = eventdata{indE(j),2}; % Event time
    ls(j,2) = find(ismember(modlist,eventdata(indE(j),14)) == 1); % Event mechanism
    ls(j,3) = aT3(indE(j));
    
    % Sets an amplitude threshold for each mechanism type
    aT = min(aT2(aT2(eventtime<mean(stress(stress(:,2)==max(stress(:,2)),1)),ls(j,2))~=0,ls(j,2)))...
        + ampthresh*(range(aT2(aT2(eventtime<mean(stress(stress(:,2)==max(stress(:,2)),1)),ls(j,2))~=0,ls(j,2))));
    
    % Seperates mechanism by P- or S-dominance
%     if ampdiff(j) < 0
%         ls2(j) = mechsep(ls(j,2));
%     elseif ampdiff(j) > 0
%         ls2(j) = mechsep(ls(j,2))+1;
%     else
%         ls2(j) = NaN;
%     end
%     
    % Seperates mechanism by amplitude
    if ls(j,3) < aT
        ls2(j) = mechsep(ls(j,2));
    else
        ls2(j) = mechsep(ls(j,2))+1;
    end 
    
    
    
end
ls(:,2) = ls2;

% Calculate probability densities
test = [];
for i = 1:max(mechsep)+1
    
    % Strain
    if TorS == 2
        gridx1 = [straintime(indE(1):nEvents:end,1);max(max(deform(:,2)))];
        x = straintime(ls(:,2) == i,1);
        try
            [f,xi,bw] = ksdensity(x,gridx1,'bandwidth',pdfsmooth*0.001);
        catch
            f = zeros(length(gridx1),1);
        end
        
    % Time
    else
        gridx1 = ls(indE(round(nEvents/2):nEvents:end),1);
        x = ls(ls(:,2) == i,1);
        try
            smoot = find(abs(deform(:,2)-(pdfsmooth*0.001)) == min(abs(deform(:,2)-(pdfsmooth*0.001))));
            [f,xi,bw] = ksdensity(x,gridx1,'bandwidth',deform(smoot(1),1));
        catch
            f = zeros(length(gridx1),1);
        end
    end
    test(:,i) = f;
end

% Converts density to a percentage
test2 = [];
for i=1:size(test,1)
    for k = 1:max(mechsep)+1
        test2(i,k) = test(i,k)/sum(test(i,:));
    end
end

% Plots percentage data
figure(10); %title(pressure);
left_color = [0 0 0]; right_color = [0 0 0];
set(gcf,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left; cla; hold on;
h = area(xi,test2.*100,'linewidth',1.1);

% Colours
h(1).FaceColor = C1;
h(2).FaceColor = brighten(C1,0.5);
h(3).FaceColor = C2;
h(4).FaceColor = brighten(C2,0.8);
h(5).FaceColor = C3;
h(6).FaceColor = brighten(C3,0.9);%

% Plot stuff
xt = get(gca, 'YTick');
set(gca, 'FontSize', 14,'Ycolor','k')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
ylabel('Relative mechanism percentage')
xlabel('Strain (%)')
set(gcf,'color','white')
ylim([0 100])

% Plots mechanical data
yyaxis right; cla; hold on;
hold on
plot(deform(:,TorS),(stress(:,2)),'-','color',[0.5 0.5 0.5],'linewidth',4);
plot(deform(:,TorS),(stress(:,2)),'k-','linewidth',3);

% Adds a box
plot([0 max(deform(:,TorS)) max(deform(:,TorS)) 0 0],...
    [0 0 1.1*max(stress(:,2)) 1.1*max(stress(:,2)) 0],'k-','linewidth',2)

% Plot stuff
ylabel('Differential stress (MPa)')
ylim([0 1.1*max(stress(:,2))])
xt = get(gca, 'YTick');
set(gca, 'FontSize', 30,'Ycolor','k')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 30)
pbaspect([4 2 1])
if TorS == 1
    xlim([0 max(deform(:,TorS))])
else
    xlim([0 max(deform(:,TorS))])
end

badj = [0.7 0.95];
B = badj.*deform(round(mean(find(max(stress(:,2))==stress(:,2)))),2);

for i = 1:length(B)
    yyaxis left
    plot([B(i) B(i)],[0 100],'w-','linewidth',4)
    plot([B(i) B(i)],[0 100],'k-','linewidth',3)
end


if saving == 1
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    %myaa('publish')
    saveas(gcf,'mechprobability.png')
end

