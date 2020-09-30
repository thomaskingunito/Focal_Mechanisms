%% Mechanism Characterisation

%% Version
% Version 1.0, 27th January 2019. Thomas King
%   - First Version

%% Parameter customisation
% Below are the suggested parameters to be modified. I don't recommend
% changing any of the code outside of these parameters.

clear all; close all

% Plotting Parameters
pressure = '5 MPa'; % This is the title of the plot
TorS = 2; % time or strain 1 or 2
averagepolarity = 0; % 1 on 0 off
ampthresh = 0.25; % Seperate mechanisms with an amplitude threshold
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
    
    % Seperates mechanism by amplitude
    if ls(j,3) < aT
        ls2(j) = mechsep(ls(j,2));
    else
        ls2(j) = mechsep(ls(j,2))+1;
    end 
end
ls(:,2) = ls2;
%%

index = 0; allsig = []; sigsamp = 4098; trainingsize = 10000;
for i = 1:length(eventdata)
        t = 0;
        while t < length(eventdata)/2
            t = randi(length(eventdata)); % Random waveform selection
        end
        try
            signal  = eventdata{i,6};
            if index == 0
                % Obtain sampling rate
                Ts = 1e-7; % Obtain sampling frequency
                Fs = 1/Ts; 
            end
        catch
            continue
        end
        [~,order] = sort(rms(signal),'descend');
        for o = 1%:size(signal,2)
            r = order(o);
            index = index+1; % Compiles data into a single matrix
            allsig(:,index) = resample(signal(:,r),sigsamp,length(signal(:,r)));
            ptimes = eventdata{i,7};
            pktimes(index) = ptimes(r);
        end
        if index > trainingsize % Stop compiling when dataset is big enough
            break
        end
        
    end
    %cd ..
    save rays_mech.mat allsig Ts % Store waveform data
    
    %%
    mechchar = [];
    for i = 1:size(allsig,2)
        csig = allsig(:,i);
        cpick = round(pktimes(i)/Ts);
        mechchar(i,1) = ls(i,2); % type
        esignal = smooth(envelope(csig,10,'rms'),25); % smoothed envelope
        try
        dsignal = inteFD(esignal(cpick:end).^2,1); % time integral
        mechchar(i,2) = (find(dsignal == max(dsignal)) - find(dsignal == min(dsignal)))*Ts; % duration
    
        % Average frequency content
        mechchar(i,3) = meanfreq(csig,Fs,[0 300000]); % whole
        
             peak = find(esignal(cpick:cpick+500) == max(esignal(cpick:cpick+500)))+cpick;
             mechchar(i,6)  = meanfreq(csig(cpick:peak),Fs,[0 300000]); % forward
     mechchar(i,7) = meanfreq(csig(peak:peak+(peak-cpick)),Fs,[0 300000]); % back
            mechchar(i,8) = obw(csig,Fs,[0 50000],25); % coda
            mechchar(i,9) = obw(csig(peak:peak+(peak-cpick)),Fs,[0 50000],25);
        
       

    [~,ang] = rangeangle([peak/Ts;esignal(peak);0],[cpick/Ts;0;0]);
    mechchar(i,4) = ang(1); % rise angle
    mechchar(i,5) = (peak-cpick)*Ts;  % peak delay
        catch
             mechchar(i,2) = NaN;
              mechchar(i,3) = NaN;
             mechchar(i,4) = NaN;
             mechchar(i,5) = NaN;
             mechchar(i,6)  = NaN;
             mechchar(i,7)  = NaN;
        end
    end
    
    