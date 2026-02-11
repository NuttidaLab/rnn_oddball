
clear;clc;
%F='D:\Neuropixels_Shah\Pixels_local\Corrected Pixels Data\';
F = '/home/shared/oddball/data/ForNuttida/Pt5';

% %Pt 4
% spikes=
% % task=
% Pt 5
%load([F,'Case05\ReanalysisAfterMC\goodData.mat']);
%load([F,'Case05\ReanalysisAfterMC\audoddball.mat']);
%load([F,'Case05\ReanalysisAfterMC\ataskpsth.mat']);
%load([F,'Case05\ReanalysisAfterMC\case05regress.mat'],'oddunits','mixedunits');

load(fullfile(F, 'goodData.mat'));
load(fullfile(F, 'audoddball.mat'));
load(fullfile(F, 'ataskpsth.mat'));
load(fullfile(F, 'case05regress.mat'), 'oddunits', 'mixedunits');

% % %Pt 6
% load([F,'Case06\CORRECTED_AP\AnalysisMC\goodData.mat']);
% load([F,'Case06\CORRECTED_AP\AnalysisMC\audodd06_ts.mat']);
% load([F,'Case06\CORRECTED_AP\AnalysisMC\ataskpsth.mat']);
% load([F,'Case06\CORRECTED_AP\AnalysisMC\case06regress.mat'],'oddunits','mixedunits');
load(fullfile(F, 'goodData.mat'));
load(fullfile(F, 'audodd06_ts.mat'));
load(fullfile(F, 'ataskpsth.mat'));
load(fullfile(F, 'case06regress.mat'), 'oddunits', 'mixedunits');

%% Set up parameters
singleU=[goodData.label]==1;
Ocells=sort([oddunits;mixedunits]);


audf1 = find(freqlabel==1);
audf2 = find(freqlabel==2);

t1std = audf1(audf1<=150);
t1odd = audf2(audf2<=150);

t2f1 = audf1(audf1>150 & audf1<=180);
t2f2 = audf2(audf2>150 & audf2<=180);

t3std = audf2(audf2>180);
t3odd = audf1(audf1>180);

o=[t1odd;t3odd]; %oddball tones
s=[t1std;t3std]; %standard tones

% ataskpsth=ataskpsth([goodData.chanDepth]<1900);
% Cell x Trial Response
nT=size(ataskpsth(1).poststim,1);
% [d,a]=deal(nan(1,nT));
for t=1:nT%is this always 330?
    for c=1:length(ataskpsth)%only oddball cells?
        %Define Response as the mean of the post stim - mean of prestim.
        %Post stim is 1.1s long, prestim is 200ms long
        R(c,t)=mean(ataskpsth(c).poststim(t,10:100)); %-mean(ataskpsth(c).prestim(t,:))*100;%10ms bins, for the 1.2 seconds after stim
    end
    p1(t)=sum(freqlabel(1:t)==1)/t;%frequenceis
    p2(t)=sum(freqlabel(1:t)==2)/t;
    n=length(max(1,t-30):t);
    p1_30(t)=sum(freqlabel(max(1,t-30):t)==1)/n;%frequenceis
    p2_30(t)=sum(freqlabel(max(1,t-30):t)==2)/n;
end

%% Distance Evolution
o1 = o(o<151);
s1 = s(s<151);
o2 = o(o>180);
s2 = s(s>180);

%[a1,d1]=deal(nan(1,length(o1)));
%X=mean(R(:,s1),2);%All standard trials
%u1=X'/norm(X);
%for t=1:length(o1)
    %u2=R(:,o1(t))';
    %u2=u2/norm(u2);
    %d1(t)=pdist2(u1,u2,"euclidean");
    %% trial(nDi,t)=t;% a(t)=pdist2(X',V(:,t)',"cosine");
    %a1(t)=Vangle(u1,u2);
%end

[a1,d1]=deal(nan(1,length(o1)));
X=mean(R(:,s2),2);%All standard trials
u1=X'/norm(X);
for t=1:length(o1)
    u2=R(:,o1(t))';
    u2=u2/norm(u2);
    d1(t)=pdist2(u1,u2,"euclidean");
    % trial(nDi,t)=t;% a(t)=pdist2(X',V(:,t)',"cosine");
    a1(t)=Vangle(u1,u2);
end

rand_trials = 1000;
d1s = zeros(rand_trials, length(o1));
for tr = 1:rand_trials
  [a1,d1]=deal(nan(1,length(o1)));
  X=mean(R(:,s2),2);%All standard trials
  u1=X'/norm(X);
  s2_sub = randsample(length(s2), length(o1), true);
  for t=1:length(s2_sub)
      u2=R(:,s2_sub(t))';
      u2=u2/norm(u2);
      d1(t)=pdist2(u1,u2,"euclidean");
      % trial(nDi,t)=t;% a(t)=pdist2(X',V(:,t)',"cosine");
      a1(t)=Vangle(u1,u2);
  end
  d1s(tr, :) = d1;
end


%[a2,d2]=deal(nan(1,length(o2)));
%X=mean(R(:,s2),2);%All standard trials
%u1=X'/norm(X);
%for t=1:length(o2)
    %u2=R(:,o2(t))';
    %u2=u2/norm(u2);
    %d2(t)=pdist2(u1,u2,"euclidean");
    %% trial(nDi,t)=t;% a(t)=pdist2(X',V(:,t)',"cosine");
    %a2(t)=Vangle(u1,u2);
%end

[a2,d2]=deal(nan(1,length(o2)));
X=mean(R(:,s1),2);%All standard trials
u1=X'/norm(X);
for t=1:length(o2)
    u2=R(:,o2(t))';
    u2=u2/norm(u2);
    d2(t)=pdist2(u1,u2,"euclidean");
    % trial(nDi,t)=t;% a(t)=pdist2(X',V(:,t)',"cosine");
    a2(t)=Vangle(u1,u2);
end

rand_trials = 1000;
d2s = zeros(rand_trials, length(o2));
for tr = 1:rand_trials
  [a2,d2]=deal(nan(1,length(o2)));
  X=mean(R(:,s1),2);%All standard trials
  u1=X'/norm(X);
  s1_sub = randsample(length(s1), length(o2), true);
  for t=1:length(s1_sub)
      u2=R(:,s1_sub(t))';
      u2=u2/norm(u2);
      d2(t)=pdist2(u1,u2,"euclidean");
      % trial(nDi,t)=t;% a(t)=pdist2(X',V(:,t)',"cosine");
      a2(t)=Vangle(u1,u2);
  end
  d2s(tr, :) = d2;
end



mdl=fitlm(1:length(d),d);
plot(mdl)
statD=[mdl.Rsquared.Ordinary,mdl.Coefficients.pValue(2)];
ylabel('Unit Vector Distance')
title(['All, r2^2=' num2str(statD(1),2) ',p=' num2str(statD(2),2)])
legend("off")
xlim([0,length(o)+1])
line(find(o<180,1,"last")*[1,1],ylim,'Color','k','LineStyle','-.')
axis tight
subplot(234);cla;
mdl=fitlm(1:length(d),a);
plot(mdl)
statD=[mdl.Rsquared.Ordinary,mdl.Coefficients.pValue(2)];
xlabel('Oddball Trial Index')
ylabel('Cosine Angle')
title(['Single Units, r2^2=' num2str(statD(1),2) ',p=' num2str(statD(2),2)])
legend("off")
xlim([0,length(o)+1])
line(find(o<180,1,"last")*[1,1],ylim,'Color','k','LineStyle','-.')
axis tight


%Single Units
subplot(232);cla;
[a,d]=deal(nan(1,length(o)));
X=mean(R(singleU,s),2);%all standard trials
u1=X'/norm(X);
for t=1:length(o)
    u2=R(singleU,o(t))';
    u2=u2/norm(u2);
    d(t)=pdist2(u1,u2,"euclidean");
    % trial(nDi,t)=t;% a(t)=pdist2(X',V(:,t)',"cosine");
    a(t)=Vangle(u1,u2);
end
mdl=fitlm(1:length(d),d);
plot(mdl)
statD=[mdl.Rsquared.Ordinary,mdl.Coefficients.pValue(2)];
xlabel('')
ylabel('')
title(['Single Units, r2^2=' num2str(statD(1),2) ',p=' num2str(statD(2),2)])
legend("off")
xlim([0,length(o)+1])
line(find(o<180,1,"last")*[1,1],ylim,'Color','k','LineStyle','-.')
axis tight
subplot(235);cla;
mdl=fitlm(1:length(d),a);
plot(mdl)
statD=[mdl.Rsquared.Ordinary,mdl.Coefficients.pValue(2)];
xlabel('Oddball Trial Index')
ylabel('')
title(['Single Units, r2^2=' num2str(statD(1),2) ',p=' num2str(statD(2),2)])
legend("off")
xlim([0,length(o)+1])
line(find(o<180,1,"last")*[1,1],ylim,'Color','k','LineStyle','-.')
axis tight

%Oddball
[a,d]=deal(nan(1,length(o)));
X=mean(R(Ocells,s),2);%all standard trials
u1=X'/norm(X);
for t=1:length(o)
        u2=R(Ocells,o(t))';
    u2=u2/norm(u2);
    d(t)=pdist2(u1,u2,"euclidean");
    % trial(nDi,t)=t;% a(t)=pdist2(X',V(:,t)',"cosine");
    a(t)=Vangle(u1,u2);
end
subplot(233);cla;
mdl=fitlm(1:length(d),d);
plot(mdl)
statD=[mdl.Rsquared.Ordinary,mdl.Coefficients.pValue(2)];
xlabel('')
ylabel('')
title(['Oddball Cells, r^2=' num2str(statD(1),2) ',p=' num2str(statD(2),2)])
legend("off")
xlim([0,length(o)+1])
line(find(o<180,1,"last")*[1,1],ylim,'Color','k','LineStyle','-.')
axis tight
subplot(236);cla;
mdl=fitlm(1:length(d),a);
plot(mdl)
statD=[mdl.Rsquared.Ordinary,mdl.Coefficients.pValue(2)];
ylabel('')
xlabel('Oddball Trial Index')
title(['Oddball Cells, r^2=' num2str(statD(1),2) ',p=' num2str(statD(2),2)])
legend("off")
xlim([0,length(o)+1])
line(find(o<180,1,"last")*[1,1],ylim,'Color','k','LineStyle','-.')
axis tight
sgtitle('Evolution of Distance, Pt 5')
%%
% %% What are the actual probabilities
% op=[p2(1:150),nan(1,30),p1(181:end)];
% op_30=[p2_30(1:150),nan(1,30),p1_30(181:end)];
% figure(1);clf;hold on;
% plot(p1);plot(p2);
% plot(p1_30);plot(p2_30)
% plot(op);plot(op_30);
% legend('Freq 1','Freq 2','Freq 1, past 30 tones','Freq 2, past 30 tones','Oddball, all tones','Oddball, past 30 tones')
% xlabel('Trial #')
% ylabel('Probability')
