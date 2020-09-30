%% Mechanism Orientations and Divergence Plots
% This code plots the orientations of the different mechanisms. The code
% stops at line 290. After, use run section for the types of plot you want.

%% Version
% Version 1.0, 27th January 2019. Thomas King
%   - First Version

%% Parameter customisation
% Below are the suggested parameters to be modified. I don't recommend
% changing any of the code outside of these parameters.

clear all; %close all; 
warning off all
% Animation and saving
anim = 0; % set to 1 to turn on animation plots for divergence maps
saving = 0; % set to 1 to save plots as they generate

% Time Corrections
timecorr = 0;

% Fracture plot options
sz = 0.4e-3; % Fracture display size

% Divergence map options
gstp = 0.005; % Gridding step

% Mechanism list
modlist = {'fitCLVD','fitDCQ','fitMM'};

% Mechanism colours
C = brighten(parula(3),.25);
C1 = C(1,:);
C2 = brighten(C(2,:),-.1);
C3 = C(3,:);

% Data windowing
wtype = 0; % Set to 0 for UCS, set to 1 for strain
badj = [0.7 0.95]; % Data windowing as a percent of UCS
wind = 0.005; % Data windowing as a value of strain
winmax = 0.1; % Maximum width of window in strain if using animation
numevents = 50; % minimum number of events per animation window

%% Compile data

% load mechanical data
stress_strain

% Load and order data
load eventdatamech_ml_residual.mat
load focalmechmodel.mat
[~,order] = sort(cell2mat(eventdata(:,2)));
eventdata = eventdata(order,:);
etime = cell2mat(eventdata(:,2)) - start - timecorr;

%% Windowing
if anim == 1
    if wtype == 1
        badj = [min(deform(:,2)):wind:max(deform(:,2))];
        B = badj;
    elseif wtype == 0
        B = badj.*deform(round(mean(find(max(stress(:,2))==stress(:,2)))),2);
    end
    for j = 1:length(B)
        B3(j) = mean(deform(abs(deform(:,2) - B(j)) == min(abs(deform(:,2) - B(j))),1));
        B(j) = find(abs(etime - B3(j)) == min(abs(etime - B3(j))));
    end
    B = unique(B);
    B = [1,B,size(eventdata,1)];
    B2 = []; index = 0; strainave = [];
    for i = 1:length(B)
        index = index+1;
        try
            stp = 1;
            B2(index,:) = [B(i-stp),B(i+stp)];
            while diff(B2(index,:)) < numevents
                B2(index,:) = [B(i-stp),B(i+stp)];
                strainave(index) = mean([badj(i-stp),badj(i+stp)]);
                stp = stp + 1;
                if diff([badj(i-stp),badj(i+stp)]) > winmax || i+stp > length(B)
                    break
                end
            end
        catch
            index = index-1;
        end
    end
    ind = find(diff(B2')>100);
    badj = badj(ind); B2 = B2(ind,:); strainave = strainave(ind(ind<length(strainave)));
    B2(1,1) = 1;
    B2(end,end) = size(eventdata,1);
else
    if wtype == 1
        badj = [0:wind:max(deform(:,2))];
        B = badj;
    elseif wtype == 0
        B = badj.*deform(round(mean(find(max(stress(:,2))==stress(:,2)))),2);
    end
    for j = 1:length(B)
        B3(j) = mean(deform(abs(deform(:,2) - B(j)) == min(abs(deform(:,2) - B(j))),1));
        B(j) = find(abs(etime - B3(j)) == min(abs(etime - B3(j))));
    end
    B = [1,B,size(eventdata,1)]; index = 0;
    for i = 1:length(B)-1
        try
            index = index + 1;
        B2(i,:) = [B(i),B(i+1)];
        strainave(i) = mean([badj(i),badj(i+1)]);
        catch
            index = index -1;
        end
    end
    ind = find(diff(B2')>100);
    B2 = B2(ind,:); strainave = strainave(ind(ind<length(strainave)));
    B2(1,1) = 1;
    B2(end,end) = size(eventdata,1);
end


%% Calculate orientations

% Orientations for each mechanism
tensile = cell(1,1,1); ind = [];
collapse = tensile;
shear = tensile;
tindex = 0; cindex = 0; sindex = 0;
tlist = ''; clist= '';
tmech = []; smech = []; cmech = [];

for i = 1:size(eventdata,1)
    
    % Skip unsolved mechanisms
    if isempty(eventdata{i,14}) == 1 || isempty(eventdata{i,2}) == 1
        continue
    end
    
    % Time correction
    eventdata{i,2} = eventdata{i,2} - start - timecorr;
    eventtime(i) = eventdata{i,2};
    
    % Fitting
    store = eventdata{i,13};
    test = cell2mat(store(5:8,:))';
    test(:,[1,3]) = 1./test(:,[1,3]);
    test2(:,1) = test(:,1).*test(:,2);
    test2(:,2) = test(:,3).*test(:,4);
    [~,order] = sort(test2(:,2),'descend');
    
    % Orientation prep
    rot = store{3,order(1)};
    rot = store{3,order(1)};
    rotangle(i) = NaN;
    azangle(i) = NaN;
    Cang = [NaN NaN];
    
    % Mechanism prep
    loc = eventdata{i,3};
    fitMOD = eval(eventdata{i,14}); % Load mechanism model
    fit = fitMOD(:,1:3);
    XYZold = fit; XYZold = XYZold'; x0=[0 0 0].'; u=[rot(1) rot(2) rot(3)].'; deg=rot(4);
    [XYZnew, R, t] = AxelRot(XYZold, deg, u, x0); fit = XYZnew';
    
    ind1 = find(fitMOD(:,3) == max(fitMOD(:,3)))';
    ind2 = find(fitMOD(:,3) == min(fitMOD(:,3)))';
    ind3 = find(fitMOD(:,2) == max(fitMOD(:,2)))';
    ind4 = find(fitMOD(:,2) == min(fitMOD(:,2)))';
    ind5 = find(fitMOD(:,1) == max(fitMOD(:,1)))';
    ind6 = find(fitMOD(:,1) == min(fitMOD(:,1)))';
    
    if mean(ismember(eventdata{i,14},modlist{3})) == 1
        etype(i) = 3;
        ind1 = []; ind2 = [];
    elseif mean(ismember(eventdata{i,14},modlist{2})) == 1
        etype(i) = 2;
        ind5 = []; ind6 = [];
    elseif mean(ismember(eventdata{i,14},modlist{1})) == 1
        etype(i) = 1;
        ind5 = []; ind6 = [];
    end
    
    % Generate surface for rotated mechanism
    p = [];
    for pp = 1:6
        try
            p = [p;(fit(eval(['ind',num2str(pp)]),:))];
        catch
            p = [p;(fit(eval(['ind',num2str(pp)]),:))'];
        end
    end
    x = p(:,1); y = p(:,2); z = p(:,3);
    xq = linspace(min(x), max (x),10);
    yq = linspace(min(y), max (y),10);
    [X,Y] = meshgrid(xq,yq);
    Z = griddata(x,y,z, X, Y, 'cubic');
    corn = [x,y,zeros(size(z,1),1)];
    
    % Generate surface for non-rotated mechanism
    p = [];
    for pp = 1:6
        try
            p = [p;fitMOD(eval(['ind',num2str(pp)]),1:3)];
        catch
            p = [p;fitMOD(eval(['ind',num2str(pp)]),1:3)'];
        end
    end
    x = p(:,1); y = p(:,2); z = p(:,3);
    xq = linspace(min(x), max (x),10);
    yq = linspace(min(y), max (y),10);
    [X1,Y1] = meshgrid(xq,yq);
    Z1 = griddata(x,y,z, X1, Y1, 'cubic');
    corn2 = [x,y,zeros(size(z,1),1)];
    
    % Calculate angular difference between the two surfaces
    [nx1 ny1 nz1] = surfnorm(X,Y,Z);
    [nx2 ny2 nz2] = surfnorm(X1,Y1,Z1+10);
    beta = acosd(dot([nx2(:),ny2(:),nz2(:)]',[nx1(:),ny1(:),nz1(:)]'));
    alpha = atan2d(nx1(:),ny1(:));
    alpha = mean(alpha(isnan(alpha)==0));
    Cang(1) = alpha;
    Cang(2) = 90-mean(beta(isnan(beta)==0));
    
    % Angular corrections
    if Cang(1) < 0
        Cang(1) = 180 + abs(diff([Cang(1) -180]));
    end
    if ismember(eventdata{i,14},modlist{3}) == 1
        Cang(1) = Cang(1) + 90; % Tensile correction
    end
    if Cang(1) > 360
        Cang(1) = 0 + abs(diff([Cang(1) 360]));
    end
    
    azangle(i) = Cang(1); % Fracture Azimuth
    rotangle(i) = 90 - abs(Cang(2)); % Fracture Dip
    
    % Set window for current event
    if badj ~= 0
        
        tind = find(B2(:,1) <= i & B2(:,2) >= i);
        
        %         tind = find(B < i);
        %         try
        %             tind = tind(end);
        %         catch
        %             tind = 1;
        %         end
    else
        tind = find(B2(:,1) <= i & B2(:,2) >= i);
    end
    
    % Generate fracture ellipsoids
    C = loc;   % center of circle
    R = 1;    % Radius of circle
    teta=0:0.01:2*pi ;
    X=R*cos(teta);
    Y=R*sin(teta) ;
    Z = zeros(size(X));
    X = X.*4;
    Y = Y./2;
    fit = [X',Y',Z'];
    
    % Rotations
    XYZold = fit; XYZold = XYZold'; x0=[0 0 0].'; u=[1 0 0].'; deg = 90;
    [XYZnew, R, t] = AxelRot(XYZold, deg, u, x0); fit = XYZnew';
    XYZold = fit; XYZold = XYZold'; x0=[0 0 0].'; u=[0 1 0].'; deg = abs(Cang(2));
    [XYZnew, R, t] = AxelRot(XYZold, deg, u, x0); fit = XYZnew';
    XYZold = fit; XYZold = XYZold'; x0=[0 0 1].'; u=[0 0 1].'; deg = 180 + azangle(i);
    [XYZnew, R, t] = AxelRot(XYZold, deg, u, x0); fit = XYZnew';
    
    % Resizing
    X = fit(:,1).*sz+loc(1); Y = fit(:,2).*sz+loc(2); Z = fit(:,3).*sz+loc(3);
    
    % Store data for individual mechanism types
    if mean(ismember(eventdata{i,14},modlist{3})) == 1
        for ttt = 1:length(tind)
            tindex = tindex+1;
            tensile{tindex,1,tind(ttt)} = X;
            tensile{tindex,2,tind(ttt)} = Y;
            tensile{tindex,3,tind(ttt)} = Z;
            tensile{tindex,4,tind(ttt)} = C3;
            tensile{tindex,5,tind(ttt)} = loc;
            tensile{tindex,6,tind(ttt)} = Cang(1);
            tensile{tindex,7,tind(ttt)} = abs(Cang(2));
        end
    elseif mean(ismember(eventdata{i,14},modlist{2})) == 1
        
        for ttt = 1:length(tind)
            sindex = sindex+1;
            shear{sindex,1,tind(ttt)} = X;
            shear{sindex,2,tind(ttt)} = Y;
            shear{sindex,3,tind(ttt)} = Z;
            shear{sindex,4,tind(ttt)} = C2;
            shear{sindex,5,tind(ttt)} = loc;
            shear{sindex,6,tind(ttt)} = Cang(1);
            shear{sindex,7,tind(ttt)} = abs(Cang(2));
        end
    elseif mean(ismember(eventdata{i,14},modlist{1})) == 1
        for ttt = 1:length(tind)
            cindex = cindex+1;
            collapse{cindex,1,tind(ttt)} = X;
            collapse{cindex,2,tind(ttt)} = Y;
            collapse{cindex,3,tind(ttt)} = Z;
            collapse{cindex,4,tind(ttt)} = C1;
            collapse{cindex,5,tind(ttt)} = loc;
            collapse{cindex,6,tind(ttt)} = Cang(1);
            collapse{cindex,7,tind(ttt)} = abs(Cang(2));
        end
    end
    
    % Removes skipped data
    ind = [ind,i];
    
    % Converts event time to strain value
    straintime(i,1) = mean(deform(abs(deform(:,1) - eventdata{i,2}) == min(abs(deform(:,1) - eventdata{i,2})),2));
    
end

eventdata = eventdata(ind,:);
straintime = straintime(ind);
etype = etype(ind);
rotangle = rotangle(ind);
azangle = azangle(ind);

%% Seperate events by azimuth according to principle shear direction

% Compile
infstore = [straintime,abs(rotangle'),etype'];
infstore2 = [straintime,azangle',etype'];

% Remove bad data
ind = find(isnan(infstore(:,2)) == 1 | isnan(infstore2(:,2)) == 1);
infstore(ind,:) = [];
infstore2(ind,:) = [];

% Shear direction
viewang = mode(round(infstore2(infstore2(:,3) == 2,2),-1)); % - 30

% Seperation and correction
infstore2(:,2) = infstore2(:,2)-viewang;
for i = 1:size(infstore2,1)
    if infstore2(i,2) > 180
        infstore2(i,2) = -180 + abs(diff([infstore2(i,2) 180]));
    elseif infstore2(i,2) < -180
        infstore2(i,2) = 180 - abs(diff([infstore2(i,2) -180]));
    end
end
ind = find(infstore2(:,2) >= -90 & infstore2(:,2) <= 90);
infstore(:,4)  = 2; % Perpendicular to shear
infstore(ind,4)  = 1; % Parallel to shear

return
%% Plot orientations vs confining pressure

figure(22);

sett = 2; conf = 40;

for f = 1:2
    
    subplot(2,2,f+sett); hold on;
    pbaspect([2 1 1])
    for a = 1:3
        e = 4-a;
        %-2+1*e
        col = eval(['C',num2str(e)]);
        data = abs(infstore(infstore(:,4) == f & infstore(:,3) == e,2));
        errorbar(conf,90-mean(data),...
            std(data(data<mean(data))),...
            std(data(data>mean(data))),'color','k','linewidth',4,...
            'marker','o','markersize',8,'markeredgecolor','k','markerfacecolor','k','capsize',8)
                errorbar(conf,90-mean(data),...
            std(data(data<mean(data)))...
            ,std(data(data>mean(data))),'color',col,'linewidth',2,...
            'marker','o','markersize',8,'markeredgecolor',col,'markerfacecolor',col)
        xlim([0 45])
        ylim([0 90])
        ylabel('Dip (\theta)')
        xlabel('Confining pressure (MPa)')
    end
end

%% Divergence Maps
figure

% Plotting stuff
windows = FindClosestFactorization(min(min([size(shear,3) size(collapse,3) size(tensile,3)])));
if windows(1) == 1 && windows(2) > 5
    windows = FindClosestFactorization(min(min([size(shear,3) size(collapse,3) size(tensile,3)]))+1);
end

for k = 2%1:min(min([size(shear,3) size(collapse,3) size(tensile,3)]))
    
    % Set current plot
    if anim == 1
        cla; hold on;
        title([num2str(mean([badj(k) badj(k+1)])),'% Strain'])
    else
        subplot(windows(1),windows(2),k);cla; hold on;
    end

    % Compile data for window
    cCol = cell2mat(collapse(:,6,k));
    cdip = cell2mat(collapse(:,7,k));
    csip  = 90-abs(cdip);
    cCol = cCol -180;
    cCol(cCol < 0) = 360 - abs(diff([cCol(cCol < 0), zeros(length(cCol(cCol < 0)),1)]'))';
    cloc = [cell2mat(tensile(:,5,k))];
    cstrike = [cell2mat(tensile(:,6,k))];
    cdip = [cell2mat(tensile(:,7,k))];
    cloc = [cloc;cell2mat(collapse(:,5,k));cell2mat(shear(:,5,k))];
    cstrike = [cstrike;cCol;cell2mat(shear(:,6,k))];
    cdip = [cdip;csip;cell2mat(shear(:,7,k))];
    
    % Removes some more bad data
    ind = find(isnan(cstrike)==1);
    cloc(ind,:) = [];
    cstrike(ind) = [];
    cdip(ind) = [];
    
    % Converts angular data into vectors
    slocs = cloc;
    posarray = []; magarray = []; magindex = 0;
    for i = 1:size(cstrike,1)
        s1  = [cstrike(i),cdip(i)];
        lr = [-gstp/2 0 0;gstp/2 0 0];
        fit = lr;
        XYZold = fit; XYZold = XYZold'; x0=[0 0 1].'; u=[0 0 1].'; deg = 90+s1(1,1)-mode(round(infstore2(infstore2(:,3) == 2,2),-1));
        [XYZnew, R, t] = AxelRot(XYZold, deg, u, x0); fit = XYZnew';
        XYZold = fit; XYZold = XYZold'; x0=[0 0 0].'; u=[0 1 0].'; deg =  s1(1,2);
        [XYZnew, R, t] = AxelRot(XYZold, deg, u, x0); fit = XYZnew';
        lr = fit+slocs(i,:);
        a = lr(lr(:,3)==max(lr(:,3)),:);
        b = slocs(i,:);
        c = lr(lr(:,3)==min(lr(:,3)),:);
        mpoint = (a+b)./2;
        magindex = magindex + 1;
        posarray(magindex,:) = [mpoint(1),mpoint(2),mpoint(3)];
        magarray(magindex,:) = [diff([a(1),c(1)]),diff([a(2),c(2)]),diff([a(3),c(3)])];
    end
    
%     %%%%%%%%%%
%     figure;
%     subplot(1,3,1) %%%%%%%%
%     
%     quiver3(posarray(:,1),posarray(:,2),posarray(:,3)...
%         ,magarray(:,1),magarray(:,2),magarray(:,3),'color','k') 
%     
%     pbaspect([4 4 10])
%     xlim([-0.02 0.02])
%     ylim([-0.02 0.02])
%     zlim([-0.05 0.05])
%    % view([viewang 0])
%     
%     %%%%%%%%%%%%
    
    % Grid 3D vector data into 2D plane
    
    gstp2 = gstp/5;
    [X3,Y3,Z3] = meshgrid(-0.02:gstp2:0.02,-0.02:gstp2:0.02,-0.05:gstp2/(5/2):0.05);
    xx = posarray(:,1);
    yy = posarray(:,2);
    zz = posarray(:,3);
    Vxx = magarray(:,1);
    Vyy = magarray(:,2);
    Vxx = zeros(length(magarray),1);
    Vzz = magarray(:,3);
    FVx = griddata(xx,yy,zz,Vxx,X3,Y3,Z3,'natural');
    FVy = griddata(xx,yy,zz,Vyy,X3,Y3,Z3,'natural');
    FVz = griddata(xx,yy,zz,Vzz,X3,Y3,Z3,'natural');
    V3x = FVx;
    V3y = FVy;
    V3z = FVz;
    
    % Calculate divergence map
    div = divergence(V3x,V3y,V3z);
    
%     %%%%%%%%%%
%     subplot(1,3,2) %%%%%%%%
%     
%     quiver3(X3,Y3,Z3...
%         ,V3x,V3y,V3z,'color','k') 
%     
%     pbaspect([4 4 10])
%     xlim([-0.02 0.02])
%     ylim([-0.02 0.02])
%     zlim([-0.05 0.02])
%     %view([viewang 0])
%     
%     %%%%%%%%%%%%
    
    % Regrids the data a bit more
    div = div(:);
    Y3 = Y3(:); Z3 = Z3(:);
    Y4 = Y3(isnan(div)==0);
    Z4 = Z3(isnan(div)==0);
    div2 = div(isnan(div)==0);
    gstp3 = gstp/50;
    [x,y] = meshgrid(-0.02:gstp3:0.02,-0.05:gstp3:0.05);
    vq = griddata(Y4,Z4,div2,x(:),y(:),'cubic');
    
    vq(vq>5e-4) = 5e-4;
    vq(vq<-5e-4) = -5e-4;
    lim = 5e-4;
    
    
    % Plotting
    contourf(-x,y,reshape(vq,size(x,1),size(x,2)),11,'LineStyle','none')
    if anim == 1
        mapstore{k,1} = strainave(k);
        mapstore{k,2} = x;
        mapstore{k,3} = y;
        mapstore{k,4} = vq;
    end
    colormap(jet);

    cmin = -lim;
    cmax = lim;
    caxis([cmin cmax])
        colorbar('location','southoutside','Ticks',[cmin,cmax],...
        'TickLabels',{'Compaction','Dilation'},'FontSize',30)
    
    % Plot bounding box
    plot([-0.02 -0.02 0.02 0.02 -0.02],[-0.05 0.05 0.05 -0.05 -0.05],'k-')
    
    % 1cm Scale bar
    plot([-0.019 -0.019 -0.009 -0.009],[-0.0479 -0.049 -0.049 -0.0479],'k-','linewidth',3)
    plot([-0.019 -0.019 -0.009 -0.009],[-0.048 -0.049 -0.049 -0.048],'w-','linewidth',1.5)
    
    % Plot stuff
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    set(gca,'box','off')
    set(gcf,'color','w')

    drawnow
    display(k)
    
    if anim == 1 && saving == 1
        cd focalmechanimation
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        %myaa(10)
        saveas(gcf,[num2str(k),'.png'])
        cd ..
    end
end

if anim == 0 && saving == 1
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    %myaa('publish')
    saveas(gcf,'divergencemaps.png')
end

%return

%% Mechanism Location plots

% Sample cylinder
[X1 Y1 Z1] = cylinder(0.02); % Makes a cylinder with radius 0.02
Z1(2,:) = 0.1; % Sets cylinder height to 0.1
shp = surf2patch(X1,Y1,Z1); % Makes it into a patch
X = shp.vertices(:,1);
Y = shp.vertices(:,2);
Z = shp.vertices(:,3)-0.05; % Puts into the correct place
shp = alphaShape(X,Y,Z,1,'HoleThreshold',10000); % Makes it into a shape

% Plot Tensile events
figure;
for k = 1:min(min([size(shear,3) size(collapse,3) size(tensile,3)]))
    
    % Plot stuff
    subplot(1,min(min([size(shear,3) size(collapse,3) size(tensile,3)])),k);hold on;
    pbaspect([4 4, 10])
    xlim([-0.02 0.02])
    ylim([-0.02 0.02])
    zlim([-0.05 0.05])
    view([viewang 0])
    
    % Plot fracture ellipses
    for i = 1:size(tensile(:,:,k),1)
        if i == 1; cla;    plot(shp,'FaceColor','black','EdgeColor','none','Facealpha',0.1)
            set(gcf,'color','w')
            axis off
            %camproj('perspective')
        end
        fill3(tensile{i,1,k},tensile{i,2,k},tensile{i,3,k},tensile{i,4,k},'linestyle','-','edgecolor',[0.25 0.25 0.25]);
    end
%     
%     % 1cm Scale bar
%     plot3([0 0 0 0],[-0.019 -0.019 -0.009 -0.009],[-0.0479 -0.049 -0.049 -0.0479],'k-','linewidth',3)
%     plot3([0 0 0 0],[-0.019 -0.019 -0.009 -0.009],[-0.048 -0.049 -0.049 -0.048],'w-','linewidth',1.5)
    
end

% Plot Shearing events
%figure;
for k = 1:min(min([size(shear,3) size(collapse,3) size(tensile,3)]))
    
    % Plot stuff
    subplot(1,min(min([size(shear,3) size(collapse,3) size(tensile,3)])),k);hold on;
    pbaspect([4 4, 10])
    xlim([-0.02 0.02])
    ylim([-0.02 0.02])
    zlim([-0.05 0.05])
    view([viewang 0])
    
    % Plot fracture ellipses
    for i = 1:size(shear(:,:,k),1)
        %         if i == 1; %cla;
        %             plot(shp,'FaceColor','black','EdgeColor','none','Facealpha',0.1)
        %             set(gcf,'color','w')
        %             axis off
        %             camproj('perspective')
        %         end
        fill3(shear{i,1,k},shear{i,2,k},shear{i,3,k},shear{i,4,k},'linestyle','-','edgecolor',[0.25 0.25 0.25]);
    end
end

% Plot Closing events
%figure(2);
for k = 1:min(min([size(shear,3) size(collapse,3) size(tensile,3)]))
    
    % Plot stuff
    subplot(1,min(min([size(shear,3) size(collapse,3) size(tensile,3)])),k);hold on;
    pbaspect([4 4, 10])
    xlim([-0.02 0.02])
    ylim([-0.02 0.02])
    zlim([-0.05 0.05])
    view([viewang 0])
    
    % Plot fracture ellipses
    for i = 1:size(collapse(:,:,k),1)
        %         if i == 1; %cla;
        %             plot(shp,'FaceColor','black','EdgeColor','none','Facealpha',0.1)
        %             set(gcf,'color','w')
        %             axis off
        %             camproj('perspective')
        %         end
        fill3(collapse{i,1,k},collapse{i,2,k},collapse{i,3,k},collapse{i,4,k},'linestyle','-','edgecolor',[0.25 0.25 0.25]);
    end
end

if saving == 1
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    %myaa('publish')
    saveas(gcf,'mechmaps.png')
end

%%


