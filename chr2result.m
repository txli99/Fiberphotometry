%CWQ,20230520
%%chr2-gcamp
clear
close all
file=dir('*.csv');
f=15;%帧率(HZ)
n=4500;%300s
k=1;
for i=1:2:length(file)
    d2=importdata(file(i).name);
    d1=importdata(file(i+1).name);
    ch1=d1.data(:,1);
    ch2=d2.data(:,1);
    
    f1=zero_filter_lowpass(ch1(1:n,1),f,2);
    f2=zero_filter_lowpass(ch2(1:n,1),f,2);
    F1=bleCorrect2([1:n]',f1,7);%多项式校正
    F2=bleCorrect2([1:n]',f2,7);%多项式校正
    total.raw_gcamp(:,k)=f1;
    total.raw_mcherry(:,k)=f2;
    ployfit_f1=F1+median(f1);
    ployfit_f2=F2+median(f2);
    total.ployfit_gcamp(:,k)=ployfit_f1;
    total.ployfit_mcherry(:,k)=ployfit_f2;
    DF_f1=(ployfit_f1-median(ployfit_f1))./median(ployfit_f1)*100;%DF/F
    DF_f2=(ployfit_f2-median(ployfit_f2))./median(ployfit_f2)*100; %DF/F
    total.df_gcamp(:,k)=DF_f1;
    total.df_mcherry(:,k)=DF_f2;
    p = polyfit(DF_f2, DF_f1, 1);
    G_fit = p(2) + p(1) *  DF_f2;
    mc_f= DF_f1- G_fit;
    total.motioncorrect(:,k)=mc_f;
    figure
    subplot(3,1,1)
    title('baseline ployfit correct')
    hold on
    plot(total.ployfit_gcamp(:,k),'g')
    plot(total.ployfit_mcherry(:,k),'m')
    subplot(3,1,2)
    title('baselien DF/F')
    hold on
    plot(total.df_gcamp(:,k),'g')
    plot(total.df_mcherry(:,k),'m')
    subplot(3,1,3)
    title('baselien motion correct DF/F')
    hold on
    plot(total.motioncorrect(:,k),'g')
    
    k=k+1;
end %导入数据
%% mean/peak df/f plot
clearvars -except total
close all
f=15;
k=1;
for j=1:size(total.motioncorrect,2)
    for i=1:1
        Light{j}(:,i)=total.motioncorrect(58*f+(i-1)*22*f+1:58*f+(i-1)*22*f+4*f,j);
        lightresponse(:,k)=Light{j}(:,i);
        k=k+1;
    end
%     meanlight(:,j)=mean(Light{j},2);
end
%mean
pre=mean(lightresponse(1:2*f,:));
light=mean(lightresponse(2*f+1:4*f,:));
post=mean(lightresponse(4*f+1:end,:));
meanresult=[pre' light' post'];
% [h,p1]=ttest(meanresult(:,1),meanresult(:,2));
figure
set(gcf,'unit','centimeters','position',[3 5 10 8 ])
set(gca,'position',[0.3,0.2,0.7,0.75] )
hold on
fillcolor=[255,255,255;147,224,255;255,255,255;]/255;
edgecolor=[147,224,255;147,224,255;147,224,255]/255;
position=[1 2 3];
m=mean(meanresult);
err=std(meanresult)./sqrt(size(meanresult,1));
for i=1:3
    bar(position(i),m(i),0.8,'FaceColor', fillcolor(i,:),'EdgeColor',edgecolor(i,:),'LineWidth',1.5)
   [x1(:,i),~]= swarmplot(position(i),meanresult(:,i),50,'k','w');
    errorbar(position(i),m(i),[],err(i),'Color',edgecolor(i,:),'Linewidth',1.5);
end
for j=[1 2]
    for ik=1:length(x1)
        line([x1(ik,j) x1(ik,j+1)],[meanresult(ik,j) meanresult(ik,j+1)],'color','k','linewidth',0.8)
    end
end
xlim([0.5,3.5]),
xticks([1:3]),xticklabels({'Pre','Light','Post'});
ylabel(['Mean' '\Delta' 'F/F(%)'])
% ylim([-0.5,1]),yticks([-0.5:0.5:3])
set(gca,'layer','top');
set(gca,'color','none','xcolor',[0 0 0],'ycolor',[0 0 0],'Linewidth',1.6,'FontName','Arial','FontSize',12,'FontWeight','bold','tickdir','out','ticklength',[0.04 0.025]);
y=ylim;
% if p1<0.05
%     ySig = max(max(meanresult(:,[1 2]))); 
%     sigline([1,2],ySig, p1(1));
% end

%peak
%mean
pre=max(lightresponse(1:2*f,:));
light=max(lightresponse(2*f+1:4*f,:));
post=max(lightresponse(4*f+1:end,:));
meanresult=[pre' light' post'];
% [h,p1]=ttest(meanresult(:,1),meanresult(:,2));
figure
set(gcf,'unit','centimeters','position',[3 4 10 8 ])
set(gca,'position',[0.3,0.2,0.7,0.75] )
hold on
fillcolor=[255,255,255;147,224,255;255,255,255;]/255;
edgecolor=[147,224,255;147,224,255;147,224,255]/255;
position=[1 2 3];
m=mean(meanresult);
err=std(meanresult)./sqrt(size(meanresult,1));
for i=1:3
    bar(position(i),m(i),0.8,'FaceColor', fillcolor(i,:),'EdgeColor',edgecolor(i,:),'LineWidth',1.5)
   [x1(:,i),~]= swarmplot(position(i),meanresult(:,i),50,'k','w');
    errorbar(position(i),m(i),[],err(i),'Color',edgecolor(i,:),'Linewidth',1.5);
end
for j=[1 2]
    for ik=1:length(x1)
        line([x1(ik,j) x1(ik,j+1)],[meanresult(ik,j) meanresult(ik,j+1)],'color','k','linewidth',0.8)
    end
end
xlim([0.5,3.5]),
xticks([1:3]),xticklabels({'Pre','Light','Post'});
ylabel(['Peak' '\Delta' 'F/F(%)'])
% ylim([-0.5,1]),yticks([-0.5:0.5:3])
set(gca,'layer','top');
set(gca,'color','none','xcolor',[0 0 0],'ycolor',[0 0 0],'Linewidth',1.6,'FontName','Arial','FontSize',12,'FontWeight','bold','tickdir','out','ticklength',[0.04 0.025]);
y=ylim;
% if p1<0.05
%     ySig = max(max(meanresult(:,[1 2]))); 
%     sigline([1,2],ySig, p1(1));
% end
%% lineplot df/f
close all
f=15;
c=[147,224,255]/255;
Lincolor=c;
data=lightresponse;
x=[-2:1/f:10];
x(1)=[];
figure
set(gcf,'unit','centimeters','position',[3 4 10 8 ])
set(gca,'position',[0.2,0.2,0.75,0.75] )
hold on
y=data;
n1=size(y,2);
m1=mean(y,2);
sd1=std(y,1,2);
sem=sd1/sqrt(n1);
shadedErrorBarsub(x,m1,sem,Lincolor(1,:));

xlabel('Time (s)')
ylabel( ['\Delta' 'F/F (%)'])

% ylim([-0.5,1]);yticks([-1:0.5:1]);
% line([0 0],[-1,1],'linestyle',':','Linewidth',1.5,'color','k')
set(gca,'color','none','xcolor',[0 0 0],'ycolor',[0 0 0],'Linewidth',1.6,'FontName','Arial','FontSize',12,'FontWeight','bold','tickdir','out')
%% 热图 df/f
close all
f=15;
c=[88,94,92;206,60,79]/255;
Lincolor=c;
pre=total.baseline.motioncorrect;
light=total.cno.motioncorrect;
h1=pre';h2=light';
if length(pre)>3000
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
