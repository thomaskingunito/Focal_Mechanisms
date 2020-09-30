

tsig_low = zeros(4096,1);
ssig_low = zeros(4096,1);
csig_low = zeros(4096,1);
tsig_high = zeros(4096,1);
ssig_high = zeros(4096,1);
csig_high = zeros(4096,1);

tsig_low_freq = zeros(1,1);
ssig_low_freq  = zeros(1,1);
csig_low_freq  = zeros(1,1);
tsig_high_freq  = zeros(1,1);
ssig_high_freq  = zeros(1,1);
csig_high_freq  = zeros(1,1);

c_amp = nan.*ones(length(ls2),5);
s_amp = nan.*ones(length(ls2),5);
t_amp = nan.*ones(length(ls2),5);
vpvs = nan.*ones(length(ls2),4);

load recloc.mat
ph = zeros(7506,1);
ampdiff = nan.*ones(length(ls2),1);
for i = 1:length(ls2)
    csig = eventdata{i,6};
%     if ismember(i,indT)
%        % continue
%     end
    
    r = rms(csig)==max(rms(csig));
    ptime = eventdata{i,7}; ptime = round(ptime(r)/1e-7);
    
    loc = eventdata{i,3};
    hyp = norm(loc-recloc(r,:));
    etime = straintime(i);   
    
    if isnan(ptime)% || etime > deform(end,1)
        continue
    end

    Vp = (hyp/(ptime*1e-7));
    
%     if Vp > 6000 < Vp < 1000
%         figure; 
%         plot(csig(:,r)); hold on;
%         scatter(ptime,0,'filled')
%         Vp
%         break
%     end
    Vs = Vp/1.3;
    stime = round((hyp/Vs)/1e-7);

    stime = 500 + (stime-ptime); 
    
    
    
    ptime = ptime-500;
    try
        csig = csig(ptime:end,r);
        if rms(csig(500:500+50))/rms(csig(1:500))<1 | max(abs(csig(500:500+100))) < 0.02 | rms(csig(500:500+50))/rms(csig(500-50:500)) > 90 
            continue
        end
    catch
        continue
    end
    
%     if i == 1
%         [csig Hd] = bandpass(csig,[1e4 6e4],1e6);
%     else
%         csig = filter(Hd,csig);
%     end
    
%     [pxx,f] = plomb(csig,1e6);
%     if ls2(i) == 6
%     ph(1:length(f)/2) = (ph(1:length(f)/2) + pxx(1:length(f)/2))/2;
%     figure(2); cla; plot(f(1:length(f)/2),ph(1:length(f)/2)); xlim([0e4 5e4])
%     title(num2str(ls2(i)))
%     drawnow
%     end
  %  if ls2(i) == 6 || ls2(i) == 2 || ls2(i) == 4
       try
        [vq1 ff] = cqt(csig(stime-200:stime+200),'BinsPerOctave',1,'SamplingFrequency',1000000,'FrequencyLimits',[3e4 6e4]);
        col = vq1.c;
%         freq = ff;
         stp = 400/size(col,2);
         col = smooth(real(col(1,:)));
         fp = findpeaks(col);
         ind = find(abs(fp-(200/stp)) == min(abs(fp-(200/stp))));
         speak = fp(ind);
         for f = 1:speak
             if col(speak-f-1) > col(speak-f)
                 break
             end
         end
         stime = (stime-200) + round(speak*stp);
%          figure(1); cla; plot(csig); hold on;
%          scatter(500,0,'filled');scatter(stime,0,'filled')
%          drawnow
   
%         [x,y] = meshgrid(1:5:3048,40000:5000:60000);
%         %[x,y] = meshgrid(1:stp:3048,ff(1:end-2));
%         %vq = col(:);
%         vq = griddata(1:stp:3048,ff(1:end-2),col,x(:),y(:),'cubic');
%         [a] = find(real(vq) > 0.001 & y(:) > 4.0e4 & y(:) < 6e4 & x(:) > 550 & x(:) < 1000);
%         x1 = x(:); y1 = y(:);
%         b = find(round(y1(a)) == max(round(y1(a))));
%         a = a(b);

   esig = envelope(csig(500:stime),100,'rms');
   ind = find(esig == max(esig))+500;
    
   esig = envelope(csig(stime:end),100,'rms');
   ind2 = find(esig == max(esig))+stime;
   grad = ind2 - ind;
   
   if grad < 30 | grad > 500
       continue
   end
   
%    esig = envelope(csig(1:end),100,'rms');
%    ind2 = find(esig(ind:end) < 0.05);
%    p = polyfit([1:ind2(1)+1].*1e-7,esig(ind:ind+ind2(1)),1);
%    grad = p(1);

   ftime = eventdata{i,7}; ftime = round(ftime(r)/1e-7);
   esig = envelope(csig,30,'rms');
   %if grad < 0
       vpvs(i,:) = [ [etime (hyp/(ftime*1e-7)) (hyp/((ftime+(stime-500))*1e-7))] grad];
       
       
       
  % else
 %      continue
  % end
       catch
           continue
       end
    if ls2(i) == 1
        %csig_low(1:4096-ptime+1) = (csig_low(1:4096-ptime+1) + csig);
        %csig_low_freq = (medfreq(csig(500:1000),1000000) + csig_low_freq)./2;
        %c_amp(i,:) = [(std(csig(500:stime-10))./sqrt(stime-10-500))  (std(csig(stime:stime+300))./sqrt(300)) rms((csig).^2) medfreq(csig(500:stime),1000000) medfreq(csig(stime:stime+300),1000000)];
    elseif ls2(i) == 2
        csig_high(1:4096-ptime+1) = (csig_high(1:4096-ptime+1) + csig);
        csig_high_freq = (medfreq(csig(500:1000),1000000) + csig_high_freq)./2;
        c_amp(i,:) = [(std(csig(500:stime-10))./sqrt(stime-10-500))  (std(csig(stime:stime+300))./sqrt(300)) rms((csig).^2) medfreq(csig(500:stime),1000000) medfreq(csig(stime:stime+300),1000000)];
    elseif ls2(i) == 3
        %ssig_low(1:4096-ptime+1) = (ssig_low(1:4096-ptime+1) + csig);
        %ssig_low_freq = (medfreq(csig(500:1000),1000000) + ssig_low_freq)./2;
        %s_amp(i,:) = [(std(csig(500:stime-10))./sqrt(stime-10-500))  (std(csig(stime:stime+300))./sqrt(300)) rms((csig).^2) medfreq(csig(500:stime),1000000) medfreq(csig(stime:stime+300),1000000)];
    elseif ls2(i) == 4
        ssig_high(1:4096-ptime+1) = (ssig_high(1:4096-ptime+1) + csig);
        ssig_high_freq = (medfreq(csig(500:1000),1000000) + ssig_high_freq)./2;
        s_amp(i,:) = [(std(csig(500:stime-10))./sqrt(stime-10-500))  (std(csig(stime:stime+300))./sqrt(300)) rms((csig).^2) medfreq(csig(500:stime),1000000) medfreq(csig(stime:stime+300),1000000)];
    elseif ls2(i) == 5
        %tsig_low(1:4096-ptime+1) = (tsig_low(1:4096-ptime+1) + csig);
        %tsig_low_freq = (medfreq(csig(500:1000),1000000) + tsig_low_freq)./2;
        %t_amp(i,:) = [(std(csig(500:stime-10))./sqrt(stime-10-500))  (std(csig(stime:stime+300))./sqrt(300)) rms((csig).^2) medfreq(csig(500:stime),1000000) medfreq(csig(stime:stime+300),1000000)];
    elseif ls2(i) == 6
        tsig_high(1:4096-ptime+1) = (tsig_high(1:4096-ptime+1) + csig);
        tsig_high_freq = (medfreq(csig(500:1000),1000000) + tsig_high_freq)./2
        t_amp(i,:) = [(std(csig(500:stime-10))./sqrt(stime-10-500))  (std(csig(stime:stime+300))./sqrt(300)) rms((csig).^2) medfreq(csig(500:stime),1000000) medfreq(csig(stime:stime+300),1000000)];
    end
    ampdiff(i) = diff([(std(csig(500:stime-10))./sqrt(stime-10-500))  (std(csig(stime:stime+300))./sqrt(300)) ]);
    
end

tsig_high = tsig_high./sum(t_amp(:,1)>0);
ssig_high = ssig_high./sum(s_amp(:,1)>0);
csig_high = csig_high./sum(c_amp(:,1)>0);

col = 3;
%%
%rotangle = 90-rotangle;
figure; hold on;
subplot(1,3,1);hold on;scatter(log(c_amp(c_amp(:,3)>0.0 & c_amp(:,1)>0.0,3)),...
    diff([(c_amp(c_amp(:,3)>0.0 & c_amp(:,1)>0.0,1)) ...
    (c_amp(c_amp(:,3)>0.0 & c_amp(:,1)>0.0,2))]'),15,real((vpvs(c_amp(:,3)>0.0 & c_amp(:,1)>0.0,4))),'o','filled');
ylim([-.05 .05])
colormap((jet(12)))
%caxis([5 12])
subplot(1,3,2);hold on;scatter(log(s_amp(s_amp(:,3)>0.0 & s_amp(:,1)>0.0,3)),...
    diff([(s_amp(s_amp(:,3)>0.0 & s_amp(:,1)>0.0,1)) ...
    (s_amp(s_amp(:,3)>0.0 & s_amp(:,1)>0.0,2))]'),15,real((vpvs(s_amp(:,3)>0.0 & s_amp(:,1)>0.0,4))),'o','filled');
ylim([-.05 .05])
%caxis([5 12])
subplot(1,3,3);hold on;scatter(log(t_amp(t_amp(:,3)>0.0 & t_amp(:,1)>0.0,3)),...
    diff([(t_amp(t_amp(:,3)>0.0 & t_amp(:,1)>0.0,1)) ...
    (t_amp(t_amp(:,3)>0.0 & t_amp(:,1)>0.0,2))]'),15,real((vpvs(t_amp(:,3)>0.0 & t_amp(:,1)>0.0,4))),'o','filled');
ylim([-.05 .05])
%caxis([5 12])
%plot(x1,y1,'r-','linewidth',2)
%%

figure; hold on;
subplot(2,3,1);hold on;scatter(log(c_amp(c_amp(:,3)>0.0 & c_amp(:,1)>0.0,3)),...
    diff([(c_amp(c_amp(:,3)>0.0 & c_amp(:,1)>0.0,1)) ...
    (c_amp(c_amp(:,3)>0.0 & c_amp(:,1)>0.0,2))]'),15,...
    vpvs(~isnan(vpvs(:,1)) & ~isnan(c_amp(:,1)),2)./vpvs(~isnan(vpvs(:,1)) & ~isnan(c_amp(:,1)),3)...
    ,'o','filled');
ylim([-.05 .05])
xlim([-6.5 0])
colormap((jet(12)))
caxis([1.2 1.4])
pbaspect([1 1 1])
subplot(2,3,2);hold on;scatter(log(s_amp(s_amp(:,3)>0.0 & s_amp(:,1)>0.0,3)),...
    diff([(s_amp(s_amp(:,3)>0.0 & s_amp(:,1)>0.0,1)) ...
    (s_amp(s_amp(:,3)>0.0 & s_amp(:,1)>0.0,2))]'),15,...
    vpvs(~isnan(vpvs(:,1)) & ~isnan(s_amp(:,1)),2)./vpvs(~isnan(vpvs(:,1)) & ~isnan(s_amp(:,1)),3)...
    ,'o','filled');
ylim([-.05 .05])
caxis([1.2 1.4])
xlim([-6.5 0])
pbaspect([1 1 1])
subplot(2,3,3);hold on;scatter(log(t_amp(t_amp(:,3)>0.0 & t_amp(:,1)>0.0,3)),...
    diff([(t_amp(t_amp(:,3)>0.0 & t_amp(:,1)>0.0,1)) ...
    (t_amp(t_amp(:,3)>0.0 & t_amp(:,1)>0.0,2))]'),15,...
    vpvs(~isnan(vpvs(:,1)) & ~isnan(t_amp(:,1)),2)./vpvs(~isnan(vpvs(:,1)) & ~isnan(t_amp(:,1)),3)...
    ,'o','filled');ylim([-.05 .05])
caxis([1.2 1.4])
ylim([-.05 .05])
xlim([-6.5 0])
pbaspect([1 1 1])

subplot(2,3,[4 5 6]);
left_color = [0 0 0]; right_color = [0 0 0];
set(gcf,'defaultAxesColorOrder',[left_color; right_color]);
hold on

plot(vpvs(~isnan(vpvs(:,1)) & ~isnan(c_amp(:,1)),1),smooth(vpvs(~isnan(vpvs(:,1))...
    & ~isnan(c_amp(:,1)),2)./vpvs(~isnan(vpvs(:,1)) & ~isnan(c_amp(:,1)),3),...
    round(sum(~isnan(c_amp(:,1)))/5)+1,'rloess')...
    ,'linewidth',4,'color','k')
plot(vpvs(~isnan(vpvs(:,1)) & ~isnan(c_amp(:,1)),1),smooth(vpvs(~isnan(vpvs(:,1))...
    & ~isnan(c_amp(:,1)),2)./vpvs(~isnan(vpvs(:,1)) & ~isnan(c_amp(:,1)),3),...
    round(sum(~isnan(c_amp(:,1)))/5)+1,'rloess')...
    ,'linewidth',1,'color',C1)
plot(vpvs(~isnan(vpvs(:,1)) & ~isnan(s_amp(:,1)),1),smooth(vpvs(~isnan(vpvs(:,1))...
    & ~isnan(s_amp(:,1)),2)./vpvs(~isnan(vpvs(:,1)) & ~isnan(s_amp(:,1)),3),...
    round(sum(~isnan(s_amp(:,1)))/5)+1,'rloess')...
    ,'linewidth',4,'color','k')
plot(vpvs(~isnan(vpvs(:,1)) & ~isnan(s_amp(:,1)),1),smooth(vpvs(~isnan(vpvs(:,1))...
    & ~isnan(s_amp(:,1)),2)./vpvs(~isnan(vpvs(:,1)) & ~isnan(s_amp(:,1)),3),...
    round(sum(~isnan(s_amp(:,1)))/5)+1,'rloess')...
    ,'linewidth',1,'color',C2)
plot(vpvs(~isnan(vpvs(:,1)) & ~isnan(t_amp(:,1)),1),smooth(vpvs(~isnan(vpvs(:,1))...
    & ~isnan(t_amp(:,1)),2)./vpvs(~isnan(vpvs(:,1)) & ~isnan(t_amp(:,1)),3),...
    round(sum(~isnan(t_amp(:,1)))/5)+1,'rloess')...
    ,'linewidth',4,'color','k')
plot(vpvs(~isnan(vpvs(:,1)) & ~isnan(t_amp(:,1)),1),smooth(vpvs(~isnan(vpvs(:,1))...
    & ~isnan(t_amp(:,1)),2)./vpvs(~isnan(vpvs(:,1)) & ~isnan(t_amp(:,1)),3),...
    round(sum(~isnan(t_amp(:,1)))/5)+1,'rloess')...
    ,'linewidth',1,'color',C3)

ylim([1.25 1.35])
xlim([min(deform(:,2)) deform(deform(:,2)==max(deform(:,2)),2)])
stress_strain
yyaxis right;
plot(deform(:,2),stress(:,2),'k-','linewidth',2)

%%
t = 1:4096;
t = t./1000;
figure;
subplot(2,3,1)
plot(t(1:length(csig_high)),csig_high,'k-')
xlim([0 1.5])

ylim([-1.1*max(abs(csig_high)) 1.1*max(abs(csig_high))])

ylabel('Amplitude')
title(['C-type'])
subplot(2,3,2)
plot(t(1:length(ssig_high)),ssig_high,'k-')
xlim([0 1.5])

ylim([-1.1*max(abs(ssig_high)) 1.1*max(abs(ssig_high))])
ylabel('Amplitude')
title(['S-type'])
subplot(2,3,3)
plot(t(1:length(tsig_high)),tsig_high,'k-')
xlim([0 1.5])
ylim([-1.1*max(abs(tsig_high)) 1.1*max(abs(tsig_high))])
title(['T-type'])

ylabel('Amplitude')


subplot(2,3,4)
%spectrogram(ssig_high(1:3048),15,[],[],1000000,'yaxis')
[vq1 ff] = cqt(csig_high(1:3048),'SamplingFrequency',1000000,'FrequencyLimits',[1e4 5e4]);
%caxis([-90 -30])
col = vq1.c;
freq = ff;
stp = 3048/size(col,2);
[x,y] = meshgrid(1:1:3048,1:100:50000);
vq = griddata(1:stp:3048,ff(1:end-2),col,x(:),y(:),'cubic');
contourf(x./1000,y./1000,reshape(real(vq),size(x,1),size(x,2)),100,'LineStyle','none')
colormap(parula)
%caxis([-0.01 0.01])
xlim([0 1.5])
ylabel('Frequency (KHz)')
xlabel('Time (ms)')
%pbaspect([1 1 1])
subplot(2,3,6)
%spectrogram(ssig_high(1:3048),15,[],[],1000000,'yaxis')
[vq1 ff] = cqt(tsig_high(1:3048),'SamplingFrequency',1000000,'FrequencyLimits',[1e4 5e4]);
%caxis([-90 -30])
col = vq1.c;
freq = ff;
stp = 3048/size(col,2);
[x,y] = meshgrid(1:1:3048,1:100:50000);
vq = griddata(1:stp:3048,ff(1:end-2),col,x(:),y(:),'cubic');
contourf(x./1000,y./1000,reshape(real(vq),size(x,1),size(x,2)),100,'LineStyle','none')
colormap(parula)
%caxis([-0.01 0.01])
xlim([0 1.5])
ylabel('Frequency (KHz)')
xlabel('Time (ms)')
%colorbar('eastoutside')
%pbaspect([1 1 1])

subplot(2,3,5)
%spectrogram(ssig_high(1:3048),15,[],[],1000000,'yaxis')
[vq1 ff] = cqt(ssig_high(1:3048),'SamplingFrequency',1000000,'FrequencyLimits',[1e4 5e4]);
%caxis([-90 -30])
col = vq1.c;
freq = ff;
stp = 3048/size(col,2);
[x,y] = meshgrid(1:1:3048,1:100:50000);
vq = griddata(1:stp:3048,ff(1:end-2),col,x(:),y(:),'cubic');
contourf(x./1000,y./1000,reshape(real(vq),size(x,1),size(x,2)),100,'LineStyle','none')
colormap(parula)
%caxis([-0.01 0.01])
xlim([0 1.5])
ylabel('Frequency (KHz)')
xlabel('Time (ms)')
%pbaspect([1 1 1])