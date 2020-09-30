%% Focal Mechanism Solutions
% Solves for focal mechanism solutions using first motion polarity
% amplitudes. Measurements are projected onto spheres and iteratively
% rotated minimising the fit to idealised mechanisms. I apologise for the
% state of 'eventdata'

clear all; close all;
warning off all

%% Initialisation

load sourceloc_ml.mat
load pktimes_ml.mat
load recloc.mat

compiledata = 0; % Set to 1 to compile waveform data
compilemodel = 0; % Set to 1 to generate fitting models

%% Compile data
if compiledata == 1
    index = 0;
    sources = cell2mat(sourcelocs(:,2:5));
    cd sg2
    for e = 1:size(sourcelocs,1)
        index = index+1;
        eventdata{index,1} = sourcelocs(e,1);
        eventdata{index,2} = sources(e,4);
        eventdata{index,3} = sources(e,1:3);
        
        clear azimuth takeoff
        for r = 1:size(recloc,1)
            
            % Azimuth and takeoff
            receiver = recloc(r,:);
            source = eventdata{index,3};
            [~,ang] = rangeangle(receiver',source');
            azimuth(r) = ang(1);
            takeoff(r) = ang(2);
            
        end
        
        eventdata{index,4} = azimuth;
        eventdata{index,5} = takeoff;
        try
            signal = leggisg2(char(eventdata{index,1}));
        catch
            index = index - 1;
            continue
        end
        eventdata{index,6} = signal;
        
        ind = find(cellfun(@isempty,pktimes(:,2,e))==1);
        pktimes(ind,2,e) = num2cell(0);
        
        eventdata{index,7} = cell2mat(pktimes(:,2,e));
        
    end
    cd ..
    save focaleventdata.mat eventdata -v7.3
else
    load focaleventdata.mat
end

%% Generate fitting spheres
if compilemodel == 1
    
    % Generate blank sphere
    [x1,y1,z1] = sphere(80);
    x = x1(:);
    y = y1(:);
    z = z1(:);
    
    % Double Couple Quad
    fitDCQ = [x y z];
    ind = find(fitDCQ(:,2) > 0 & fitDCQ(:,3) > 0);
    fitDCQ(ind,4) = 1;
    ind = find(fitDCQ(:,2) < 0 & fitDCQ(:,3) < 0);
    fitDCQ(ind,4) = 1;
    ind = find(fitDCQ(:,2) < 0 & fitDCQ(:,3) > 0);
    fitDCQ(ind,4) = -1;
    ind = find(fitDCQ(:,2) > 0 & fitDCQ(:,3) < 0);
    fitDCQ(ind,4) = -1;
    
    for i = 1:size(fitDCQ,1)
        [Idx,D] = knnsearch(fitDCQ(i,1:3),fitDCQ(:,1:3));
        ind = find(D < 0.5);
        val = fitDCQ(ind,4);
        val = val(find(isnan(val) == 0));
        fitDCQ(i,5) = mean(val);
    end
    
    figure(1); subplot_tight(1,3,2,[0.04,0.01]); surf(x1,y1,z1,reshape(fitDCQ(:,5),size(x1,1),size(x1,2)),'linestyle','none');
    set(gcf,'color','white')
    pbaspect([1 1 1])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    set(gca,'box','off')
    set(gca,'ztick',[])
    set(gca,'zticklabel',[])
    title('S-type')
    axis off
    
    % CLVD
    fitCLVD = [x y z];
    [X1 Y1 Z1] = cylinder(0.75);
    Z1(2,:) = 1;
    shp = surf2patch(X1,Y1,Z1);
    X1 = shp.vertices(:,1);
    Y1 = shp.vertices(:,2);
    Z1 = shp.vertices(:,3)-1;
    ind = find(Z1 == 0);
    X1 = vertcat(X1,X1(ind));
    Y1 = vertcat(Y1,Y1(ind));
    Z1 = vertcat(Z1,Z1(ind)+1);
    XYZold = [X1 Y1 Z1]; XYZold = XYZold'; x0=[0 0 0].'; u=[0 1 0].'; deg=90;
    [XYZnew, R, t] = AxelRot(XYZold, deg, u, x0); XYZnew = XYZnew';
    X1 = XYZnew(:,1); Y1 = XYZnew(:,2); Z1 = XYZnew(:,3);
    shp = alphaShape(X1,Y1,Z1,1,'HoleThreshold',15);
    ind = find(inShape(shp,fitCLVD(:,1),fitCLVD(:,2),fitCLVD(:,3)) == 1);
    fitCLVD(:,4) = 1;
    fitCLVD(ind,4) = -1;
    
    for i = 1:size(fitCLVD,1)
        [Idx,D] = knnsearch(fitCLVD(i,1:3),fitCLVD(:,1:3));
        ind = find(D < 0.5);
        val = fitCLVD(ind,4);
        val = val(find(isnan(val) == 0));
        fitCLVD(i,5) = mean(val);
    end
    
    figure(1); subplot_tight(1,3,1,[0.04,0.01]); surf(x1,y1,z1,reshape(fitCLVD(:,5),size(x1,1),size(x1,2)),'linestyle','none');
    view([0 30]); pbaspect([1 1 1])
    set(gcf,'color','white')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    set(gca,'box','off')
    set(gca,'ztick',[])
    set(gca,'zticklabel',[])
    title('C-type')
    axis off
    
    % Mixed Mode
    fitMM = [x y z];
    ind = find(fitMM(:,3) < 0.1 & fitMM(:,3) > -0.1);
    fitMM(:,4) = 1;
    fitMM(ind,4) = -1;
    
    for i = 1:size(fitMM,1)
        [Idx,D] = knnsearch(fitMM(i,1:3),fitMM(:,1:3));
        ind = find(D < 0.5);
        val = fitMM(ind,4);
        val = val(find(isnan(val) == 0));
        fitMM(i,5) = mean(val);
    end
    
    figure(1); subplot_tight(1,3,3,[0.04,0.01]); surf(x1,y1,z1,reshape(fitMM(:,5),size(x1,1),size(x1,2)),'linestyle','none');
    view([0 30]); pbaspect([1 1 1])
    title('T-type')
    set(gcf,'color','white')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    set(gca,'box','off')
    set(gca,'ztick',[])
    set(gca,'zticklabel',[])
    axis off
    
    save focalmechmodel.mat fitDCQ fitCLVD fitMM
end

%% Data prep
clear all; close all

load focaleventdata.mat
load focalmechmodel.mat

Fs = 1000000;           % Sampling frequency
T = 0.1/(Fs);           % Sampling period
L = 2048;               % Length of signal
rayt = (0:L-1)*T;       % Time vector

Wn = ([0.0001 20]/1000); %frequency band
[z,p,k] = butter(6,Wn,'stop'); %butter filter
[sos,g] = zp2sos(z,p,k);      % Convert to SOS form
Hd3 = dfilt.df2tsos(sos,g);   % Create a dfilt object

%% Pick first motions

fmstore = [];
for e = 1:size(eventdata,1)
    display(num2str(size(eventdata,1) - e +1))
    signal = eventdata{e,6};
    ptimes = eventdata{e,7};
    if length(ptimes) == size(signal,2)
        [~,order] = sort(max(abs(signal)),'descend');
        amp = zeros(1,length(order)); error = amp; pol = amp;
        for o = 1:size(signal,2)
            r = order(o);
            csignal = filter(Hd3,signal(:,r)); csignal = csignal - mean(csignal);
            cpick = round((ptimes(r))/T); %*10
            if cpick-100 <= 0 || isnan(cpick) == 1
                amp(r) = 0;
                pol(r) = 0;
                error(r) = 0;
                continue
            end
            env = smooth(envelope(csignal,1,'rms'),5);
            noise = max(env(1:cpick));
            cropenv = env(cpick-100:cpick+100);
            scsig = smooth(csignal,5);
            for i = 1:length(csignal)-cpick
                if i < cpick
                    continue
                end
                if env(i) > 1.1*noise && env(i+3) < env(i)
                    amp(r) = csignal(i);
                    if amp(r) < 0; pol(r) = -1; else; pol(r) = 1; end
                    error(r) = abs(diff([i cpick]));
                    break
                end
            end
        end
        eventdata{e,8} = amp; % First motion amplitude
        eventdata{e,9} = pol; % First motion polarity
        eventdata{e,10} = error; % Difference of onset to first motion pick
        fmstore = [fmstore,amp];
    end
end


%% Compile polarity data

pol = [];
for e= 1:size(eventdata,1)
    cpol = eventdata{e,9};
    if isempty(cpol) == 0
        pol = [pol;cpol];
    else
        pol = [pol;zeros(1,length(order)).*NaN];
    end
end

% Reset
sources = cell2mat(eventdata(:,3));
chk = [];
e = 1;

%% Mechanism Inversion
tic
while e ~= size(eventdata,1)+1
    display(num2str(size(eventdata,1) - e +1))
    
    % Minimum of 8 polarity measurements
    if sum(abs(pol(e,pol(e,:) ~= 0))) < 8
        e = e+1;
        continue
    end
    
    ind = e;
    eventdata{e,11} = ind;
    
    % Compile measurements
    th = cell2mat(eventdata(ind,4));
    rho = cell2mat(eventdata(ind,5));
    amp = cell2mat(eventdata(ind,8));
    th = reshape(th,size(th,1)*size(th,2),1);
    rho = reshape(rho,size(rho,1)*size(rho,2),1);
    amp = reshape(amp,size(amp,1)*size(amp,2),1);
    amp = amp./max(abs(amp));
    
    % Models to test
    global modlist
    modlist = {'fitCLVD','fitMM','fitDCQ'};
    
    % Compile input data
    [x,y,z] = sph2cart(th,rho,1);
    modS = [x(:) y(:) z(:) amp zeros(length(x(:)),1)]; % Measured
    global models
    for x = 1:size(modlist,2)
        models{x} = eval(modlist{x}); % Modeled
    end
    
    % Reset
    storebck = cell(11,size(modlist,2));
    residualbck = [];
    mechsolbck = [];
    
    % Inversion
    for x = 1:length(modlist)
        [output] = focmech(modS,x); % x is current model
        store = output{3};
        storebck(:,x) = store(:,x); % Inversion results
        residualbck(x) = output{2}; % Residual
        
    end
    store = storebck;
    eventdata{e,13} = store;
    
    % Mechanism is model with lowest residual
    try
        order = find(residualbck == min(residualbck));
        
        eventdata{e,14} = store{1,order(1)};
    catch
        eventdata{e,14} = [];
    end
    
    display(eventdata{e,14})
    e = e+1;
    
    % Autosave
    if toc > 3000
        save eventdatamech_ml_residual.mat eventdata -v7.3
        tic
    end
    
    % Event percentages
    figure(12); cla;
    A = eventdata(:,14);
    ind = find(~cellfun(@isempty,A));
    A = A(ind);
    [u,~,n] = unique(A(:));
    B = accumarray(n, 1, [], @sum);
    bar(B)
    set(gca,'XTickLabel',u)
    drawnow
end

save eventdatamech_ml_residual.mat eventdata -v7.3

