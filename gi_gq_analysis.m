%CWQ,20230520
%%GI/GQ
clear
close all
file1=dir('*pre*.csv');
file2=dir('*post*.csv');
f=15;%帧率(HZ)
n=4500;%300s
k=1;
for i=1:2:length(file1)
    d2=importdata(file1(i).name);
    d1=importdata(file1(i+1).name);
    d4=importdata(file2(i).name);
    d3=importdata(file2(i+1).name);
    baseline.ch1=d1.data(:,1);
    baseline.ch2=d2.data(:,1);
    cno.ch1=d3.data(:,1);
    cno.ch2=d4.data(:,1);
    
    f1=zero_filter_lowpass(baseline.ch1(1:n,1),f,2);
    f2=zero_filter_lowpass(baseline.ch2(1:n,1),f,2);
    F1=bleCorrect2([1:n]',f1,4);%多项式校正
    F2=bleCorrect2([1:n]',f2,4);%多项式校正
    total.baseline.raw_gcamp(:,k)=f1;
    total.baseline.raw_mcherry(:,k)=f2;
    ployfit_f1=F1+median(f1);
    ployfit_f2=F2+median(f2);
    total.baseline.ployfit_gcamp(:,k)=ployfit_f1;
    total.baseline.ployfit_mcherry(:,k)=ployfit_f2;
    DF_f1=(ployfit_f1-median(ployfit_f1))./median(ployfit_f1)*100;%DF/F
    DF_f2=(ployfit_f2-median(ployfit_f2))./median(ployfit_f2)*100; %DF/F
    total.baseline.df_gcamp(:,k)=DF_f1;
    total.baseline.df_mcherry(:,k)=DF_f2;
    p = polyfit(DF_f2, DF_f1, 1);
    G_fit = p(2) + p(1) *  DF_f2;
    mc_f= DF_f1- G_fit;
    total.baseline.motioncorrect(:,k)=mc_f;
    figure
    subplot(3,2,1)
    title('baseline ployfit correct')
    hold on
    plot(total.baseline.ployfit_gcamp(:,k),'g')
    plot(total.baseline.ployfit_mcherry(:,k),'m')
    subplot(3,2,3)
    title('baselien DF/F')
    hold on
    plot(total.baseline.df_gcamp(:,k),'g')
    plot(total.baseline.df_mcherry(:,k),'m')
    subplot(3,2,5)
    title('baselien motion correct DF/F')
    hold on
    plot(total.baseline.motioncorrect(:,k),'g')
    
    f1=zero_filter_lowpass(cno.ch1(1:n,1),f,2);
    f2=zero_filter_lowpass(cno.ch2(1:n,1),f,2);
    F1=bleCorrect([1:n]',f1,4);%多项式校正
    F2=bleCorrect([1:n]',f2,4);%多项式校正
    total.cno.raw_gcamp(:,k)=f1;
    total.cno.raw_mcherry(:,k)=f2;
    ployfit_f1=F1+median(f1);
    ployfit_f2=F2+median(f2);
    total.cno.ployfit_gcamp(:,k)=ployfit_f1;
    total.cno.ployfit_mcherry(:,k)=ployfit_f2;
    DF_f1=(ployfit_f1-median(ployfit_f1))./median(ployfit_f1)*100;%DF/F
    DF_f2=(ployfit_f2-median(ployfit_f2))./median(ployfit_f2)*100; %DF/F
    total.cno.df_gcamp(:,k)=DF_f1;
    total.cno.df_mcherry(:,k)=DF_f2;
    p = polyfit(DF_f2, DF_f1, 1);
    G_fit = p(2) + p(1) *   DF_f2;
    mc_f= DF_f1- G_fit;
    total.cno.motioncorrect(:,k)=mc_f;
    
    subplot(3,2,2)
    title('cno ployfit correct')
    hold on
    plot(total.cno.ployfit_gcamp(:,k),'g')
    plot(total.cno.ployfit_mcherry(:,k),'m')
    subplot(3,2,4)
    title('cno DF/F')
    hold on
    plot(total.cno.df_gcamp(:,k),'g')
    plot(total.cno.df_mcherry(:,k),'m')
    subplot(3,2,6)
    title('cno motion correct DF/F')
    hold on
    plot(total.cno.motioncorrect(:,k),'g')
    k=k+1;
end %导入数据
%% mean/peak df/f plot
clearvars -except total
close all
f=15;
%mean
baseline=mean(total.baseline.motioncorrect);
cno=mean(total.cno.motioncorrect);
meanresult=[baseline' cno'];
[h,p1]=ttest(meanresult(:,1),meanresult(:,2));
figure
set(gcf,'unit','centimeters','position',[10 10 6 8 ])
set(gca,'position',[0.3,0.2,0.7,0.75] )
hold on
fillcolor=[88,94,92;206,60,79]/255;
edgecolor=[88,94,92;206,60,79]/255;
position=[1 2 ];
m=mean(meanresult);
err=std(meanresult)./sqrt(size(meanresult,1));
for i=1:2
    bar(position(i),m(i),0.8,'FaceColor', fillcolor(i,:),'EdgeColor',edgecolor(i,:),'LineWidth',1.5)
   [x1(:,i),~]= swarmplot(position(i),meanresult(:,i),50,'k','w');
    errorbar(position(i),m(i),[],err(i),'Color',edgecolor(i,:),'Linewidth',1.5);
end

for ik=1:length(x1)
    line([x1(ik,1) x1(ik,2)],[meanresult(ik,1) meanresult(ik,2)],'color','k','linewidth',0.8)
end

xlim([0.5,2.5]),
xticks([1 2]),xticklabels({'Baseline','CNO'});
ylabel(['Mean' '\Delta' 'F/F(%)'])
% ylim([-0.5,1]),yticks([-0.5:0.5:3])
set(gca,'layer','top');
set(gca,'color','none','xcolor',[0 0 0],'ycolor',[0 0 0],'Linewidth',1.6,'FontName','Arial','FontSize',12,'FontWeight','bold','tickdir','out','ticklength',[0.04 0.025]);
y=ylim;
if p1<0.05
    ySig = max(max(meanresult(:,[1 2]))); 
    sigline([1,2],ySig, p1(1));
end

%peak
baseline=max(total.baseline.motioncorrect);
cno=max(total.cno.motioncorrect);
peakresult=[baseline' cno'];
[h,p2]=ttest(peakresult(:,1),peakresult(:,2));
figure
set(gcf,'unit','centimeters','position',[10 10 6 8 ])
set(gca,'position',[0.3,0.2,0.7,0.75] )
hold on

position=[1 2 ];
m=mean(peakresult);
err=std(peakresult)./sqrt(size(peakresult,1));
for i=1:2
    bar(position(i),m(i),0.8,'FaceColor', fillcolor(i,:),'EdgeColor',edgecolor(i,:),'LineWidth',1.5)
   [x1(:,i),~]= swarmplot(position(i),peakresult(:,i),50,'k','w');
    errorbar(position(i),m(i),[],err(i),'Color',edgecolor(i,:),'Linewidth',1.5);
end

for ik=1:length(x1)
    line([x1(ik,1) x1(ik,2)],[peakresult(ik,1) peakresult(ik,2)],'color','k','linewidth',0.8)
end

xlim([0.5,2.5]),
xticks([1 2]),xticklabels({'Baseline','CNO'});
ylabel(['Peak' '\Delta' 'F/F(%)'])
% ylim([-0.5,1]),yticks([-0.5:0.5:3])
set(gca,'layer','top');
set(gca,'color','none','xcolor',[0 0 0],'ycolor',[0 0 0],'Linewidth',1.6,'FontName','Arial','FontSize',12,'FontWeight','bold','tickdir','out','ticklength',[0.04 0.025]);
y=ylim;
if p2<0.05
    ySig = max(max(peakresult(:,[1 2]))); 
    sigline([1,2],ySig, p2(1));
end
%% lineplot df/f
f=15;
c=[88,94,92;206,60,79]/255;
Lincolor=c;
baseline=total.baseline.motioncorrect;
cno=total.cno.motioncorrect;
x=[0:1/(f*60):5];
x(1)=[];
figure
set(gcf,'unit','centimeters','position',[10 10 30 8 ])
set(gca,'position',[0.2,0.2,0.75,0.75] )
hold on
y=baseline;
n1=size(y,2);
m1=mean(y,2);
sd1=std(y,1,2);
sem=sd1/sqrt(n1);
shadedErrorBarsub(x,m1,sem,Lincolor(1,:));

y=cno;
n2=size(y,1);
m2=mean(y,2);
sd2=std(y,1,2);
sem=sd2/sqrt(n2);
shadedErrorBarsub(x,m2,sem,Lincolor(2,:));
xlim([0,5])
xticks([0:1:5]);
xlabel('Time (min)')
ylabel( ['\Delta' 'F/F (%)'])

% ylim([-0.5,1]);yticks([-1:0.5:1]);
% line([0 0],[-1,1],'linestyle',':','Linewidth',1.5,'color','k')
set(gca,'color','none','xcolor',[0 0 0],'ycolor',[0 0 0],'Linewidth',1.6,'FontName','Arial','FontSize',12,'FontWeight','bold','tickdir','out')
%% 热图 df/f
f=15;
c=[88,94,92;206,60,79]/255;
Lincolor=c;
baseline=total.baseline.motioncorrect;
cno=total.cno.motioncorrect;

baseline=total.baseline.motioncorrect;
cno=total.cno.motioncorrect;
h1=baseline';h2=cno';
if length(baseline)>3000
    n=2*f;%%downsample
    H1=[];H2=[];
    rowD=ones(1,length(h1)/n)*n;%downsample,n row mean
    temp=mat2cell(h1',rowD);
    dstemp=[];
    for ik=1:length(temp)
        dstemp(ik,:)=mean(temp{ik});%downsample
    end
    H1=[H1 dstemp'];
    temp=mat2cell(h2',rowD);
    dstemp=[];
    for ik=1:length(temp)
        dstemp(ik,:)=mean(temp{ik});%downsample
    end
    H2=[H2 dstemp'];
    H1=movmean(H1,10,2);
    H2=movmean(H2,10,2);
end

figure
set(gcf,'unit','centimeters','position',[10 10 20 8 ])
set(gca,'position',[0.1,0.2,0.7,0.7] )
H=heatmap(H1);
H.Colormap=parula;
% H.ColorScaling='scaledrows'
H.GridVisible='off'
H.ColorLimits = [0 0.5];
for j=1:length(H.XDisplayLabels)
    H.XDisplayLabels{j}=' ';
end
% for j=1:length(H.YDisplayLabels)
%     H.YDisplayLabels{j}=' ';
% end
H.ColorbarVisible ='off';
H.MissingDataColor = [1 1 1];
figure
set(gcf,'unit','centimeters','position',[10 10 20 8 ])
set(gca,'position',[0.1,0.2,0.7,0.7] )
H=heatmap(H2);
H.Colormap=parula;
% H.ColorScaling='scaledrows';%min=0,max=1
H.GridVisible='off';
% H.ColorbarVisible ='off';
H.ColorLimits = [0 0.5];
for j=1:length(H.XDisplayLabels)
    H.XDisplayLabels{j}=' ';
end
% for j=1:length(H.YDisplayLabels)
%     H.YDisplayLabels{j}=' ';
% end
H.MissingDataColor = [1 1 1];
H.GridVisible='off';
% H.ColorbarVisible ='off';
