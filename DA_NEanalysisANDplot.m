%% total calculate,plotfigure
%%total
clear
close all
%%DA:GREEN,NE:BLUE
original=dir('*original*.csv');
control=dir('*control*.csv');
event=dir('*events*.csv');
group1=[3,4,5,28,29,30];%1DA,2NE
group2=[6:28];%1NE,2DA
base=2;
for i=1:length(event)
    D1=importdata(original(i).name);
    D2=importdata(control(i).name);
    D3=importdata(event(i).name);
    timestrp=D1.data(:,3);%time
    f=round(1/(timestrp(2,1)-timestrp(1,1))*10^9);%frame rate
    n=length(D1.data(:,1));
    timestrp=timestrp(1:n-1,:);%time
    miceid=str2num(original(i).name(1:2));
    if ~isnan(intersect(miceid,group1))
        gcamp.raw_DA{i}= medfilt1(D1.data(1:n-1,1),3);
        mcherry.raw_DA{i}= medfilt1(D2.data(1:n-1,1),3);
        gcamp.raw_NE{i}= medfilt1(D1.data(1:n-1,2),3);
        mcherry.raw_NE{i}= medfilt1(D2.data(1:n-1,2),3);
    else
        gcamp.raw_DA{i}=medfilt1(D1.data(1:n-1,2),3);
        mcherry.raw_DA{i}=medfilt1(D2.data(1:n-1,2),3);
        gcamp.raw_NE{i}=medfilt1(D1.data(1:n-1,1),3);
        mcherry.raw_NE{i}=medfilt1(D2.data(1:n-1,1),3);
    end
    % filter and polyfit   
    F1=zero_filter_lowpass(gcamp.raw_DA{i},f,2);
    F2=zero_filter_lowpass(mcherry.raw_DA{i},f,2);
    F3=zero_filter_lowpass(gcamp.raw_NE{i},f,2);
    F4=zero_filter_lowpass(mcherry.raw_NE{i},f,2);
    figure
    set(gcf,'unit','centimeters','position',[0 0 40 20])%7.25等同 graphpad 5cm
    subplot(4,1,1)
    title('raw data')
    hold on
    plot(F1,'g');plot(F2,'m'); plot(F3,'b'); plot(F4,'k');
    
    cf1=bleCorrect2(timestrp,F1,10)+median(F1(1:20*f));%多项式校正
    cf2=bleCorrect2(timestrp,F2,10)+median(F2(1:20*f));%多项式校正
    cf3=bleCorrect2(timestrp,F3,10)+median(F3(1:20*f));%多项式校正
    cf4=bleCorrect2(timestrp,F4,10)+median(F4(1:20*f));%多项式校正
    subplot(4,1,2)
    title('correct data')
    hold on
    plot(cf1,'g');plot(cf2,'m'); plot(cf3,'b'); plot(cf4,'k');
    gcamp.correct_DA{i}=cf1;
    mcherry.correct_DA{i}=cf2;
    gcamp.correct_NE{i}=cf3;
    mcherry.correct_NE{i}=cf4;
    markertime=D3.data;
    temp1=find(timestrp<=markertime(1));
    eventpoint1=temp1(end);
    temp2=find(timestrp>=markertime(2));
    eventpoint2=temp2(1);
    baseline1g=gcamp.correct_DA{i}(eventpoint1-base*60*f:eventpoint1-1,:);%5 min baseline
    baseline2g=gcamp.correct_NE{i}(eventpoint1-base*60*f:eventpoint1-1,:);%5 min baseline
    baseline1m=mcherry.correct_DA{i}(eventpoint1-base*60*f:eventpoint1-1,:);%5 min baseline
    baseline2m=mcherry.correct_NE{i}(eventpoint1-base*60*f:eventpoint1-1,:);%5 min baseline
    constranit1g=gcamp.correct_DA{i}(eventpoint2:eventpoint2+5.9*60*60*f,:);%6h
    constranit2g=gcamp.correct_NE{i}(eventpoint2:eventpoint2+5.9*60*60*f,:);%6h
    constranit1m=mcherry.correct_DA{i}(eventpoint2:eventpoint2+5.9*60*60*f,:);%6h
    constranit2m=mcherry.correct_NE{i}(eventpoint2:eventpoint2+5.9*60*60*f,:);%6h
    base_DF1=(baseline1g-median(baseline1g))./median(baseline1g)*100;
    base_DF2=(baseline2g-median(baseline2g))./median(baseline2g)*100;
    base_DF3=(baseline1m-median(baseline1m))./median(baseline1m)*100;
    base_DF4=(baseline2m-median(baseline2m))./median(baseline2m)*100;
    constr_DF1=(constranit1g-median(baseline1g))./median(baseline1g)*100;
    constr_DF2=(constranit2g-median(baseline2g))./median(baseline2g)*100;
    constr_DF3=(constranit1m-median(baseline1m))./median(baseline1m)*100;
    constr_DF4=(constranit2m-median(baseline2m))./median(baseline2m)*100;
    gcamp.df_DA{i}=[base_DF1;constr_DF1];
    gcamp.df_NE{i}=[base_DF2;constr_DF2];
    mcherry.df_DA{i}=[base_DF3;constr_DF3];
    mcherry.df_NE{i}=[base_DF4;constr_DF4];
    subplot(4,1,3)
    title('DF/F') 
    hold on
    plot(gcamp.df_DA{i},'g');plot(mcherry.df_DA{i},'m'); plot(gcamp.df_NE{i},'b'); plot(mcherry.df_NE{i},'k');
    y=ylim;
    rectangle('Position',[0,y(1),5*60*f,y(2)-y(1)],'EdgeColor','r','FaceColor',[0.8,0.8,0.8])
    hold on
    plot(gcamp.df_DA{i},'g');plot(mcherry.df_DA{i},'m'); plot(gcamp.df_NE{i},'b'); plot(mcherry.df_NE{i},'k');
    % motion correct
    p = polyfit(mcherry.df_DA{i},gcamp.df_DA{i}, 1);
    G_fit = p(2) + p(1) *  mcherry.df_DA{i};
    F1=gcamp.df_DA{i}- G_fit;
    p = polyfit(mcherry.df_NE{i},gcamp.df_NE{i}, 1);
    G_fit = p(2) + p(1) *  mcherry.df_DA{i};
    F2=gcamp.df_NE{i}- G_fit;
    subplot(4,1,4)
    title('motion correct')
    hold on
    plot(F1,'g');plot(F2,'b');
    y=ylim;
    rectangle('Position',[0,y(1),5*60*f,y(2)-y(1)],'EdgeColor','r','FaceColor',[0.8,0.8,0.8])
    hold on
    plot(F1,'g');plot(F2,'b');
    
    gcamp.df_mc_DA{i}=F1;
    gcamp.df_mc_NE{i}=F2;
    total.DF_DA_g(:,i)=gcamp.df_DA{i};
    total.DF_NE_g(:,i)=gcamp.df_NE{i};
    total.DF_DA_m(:,i)=mcherry.df_DA{i};
    total.DF_NE_m(:,i)=mcherry.df_NE{i};
    total.correctDF_DA(:,i)=gcamp.df_mc_DA{i};
    total.correctDF_NE(:,i)=gcamp.df_mc_NE{i};
    name=original(i).name(1:3);
    total.miceID(i,:)=str2num(name);
    suptitle(name)
    saveas(gcf,name,'png')
    clearvars -except total gcamp  i mcherry original control event group1 group2 base

end
%% total calculate,no figure
clear
close all
%%DA:GREEN,NE:BLUE
original=dir('*original*.csv');
control=dir('*control*.csv');
event=dir('*events*.csv');
group1=[3,4,5,28,29,30];%1DA,2NE
group2=[6:28];%1NE,2DA
base=2; %baseline length
for i=1:length(event)
    D1=importdata(original(i).name);
    D2=importdata(control(i).name);
    D3=importdata(event(i).name);
    timestrp=D1.data(:,3);%time
    f=round(1/(timestrp(2,1)-timestrp(1,1))*10^9);%frame rate
    n=length(D1.data(:,1));
    timestrp=timestrp(1:n-1,:);%time
    miceid=str2num(original(i).name(1:2));
    if ~isnan(intersect(miceid,group1))
        gcamp.raw_DA{i}= medfilt1(D1.data(1:n-1,1),3);
        mcherry.raw_DA{i}= medfilt1(D2.data(1:n-1,1),3);
        gcamp.raw_NE{i}= medfilt1(D1.data(1:n-1,2),3);
        mcherry.raw_NE{i}= medfilt1(D2.data(1:n-1,2),3);
    else
        gcamp.raw_DA{i}=medfilt1(D1.data(1:n-1,2),3);
        mcherry.raw_DA{i}=medfilt1(D2.data(1:n-1,2),3);
        gcamp.raw_NE{i}=medfilt1(D1.data(1:n-1,1),3);
        mcherry.raw_NE{i}=medfilt1(D2.data(1:n-1,1),3);
    end
    % filter and polyfit   
    F1=zero_filter_lowpass(gcamp.raw_DA{i},f,2);
    F2=zero_filter_lowpass(mcherry.raw_DA{i},f,2);
    F3=zero_filter_lowpass(gcamp.raw_NE{i},f,2);
    F4=zero_filter_lowpass(mcherry.raw_NE{i},f,2);
    
    
    cf1=bleCorrect2(timestrp,F1,10)+median(F1);%多项式校正
    cf2=bleCorrect2(timestrp,F2,10)+median(F2);%多项式校正
    cf3=bleCorrect2(timestrp,F3,10)+median(F3);%多项式校正
    cf4=bleCorrect2(timestrp,F4,10)+median(F4);%多项式校正
    
    gcamp.correct_DA{i}=cf1;
    mcherry.correct_DA{i}=cf2;
    gcamp.correct_NE{i}=cf3;
    mcherry.correct_NE{i}=cf4;
    markertime=D3.data;
    temp1=find(timestrp<=markertime(1));
    eventpoint1=temp1(end);
    temp2=find(timestrp>=markertime(2));
    eventpoint2=temp2(1);
    baseline1g=gcamp.correct_DA{i}(eventpoint1-base*60*f:eventpoint1-1,:);%5 min baseline
    baseline2g=gcamp.correct_NE{i}(eventpoint1-base*60*f:eventpoint1-1,:);%5 min baseline
    baseline1m=mcherry.correct_DA{i}(eventpoint1-base*60*f:eventpoint1-1,:);%5 min baseline
    baseline2m=mcherry.correct_NE{i}(eventpoint1-base*60*f:eventpoint1-1,:);%5 min baseline
    constranit1g=gcamp.correct_DA{i}(eventpoint2:eventpoint2+5.9*60*60*f,:);%6h
    constranit2g=gcamp.correct_NE{i}(eventpoint2:eventpoint2+5.9*60*60*f,:);%6h
    constranit1m=mcherry.correct_DA{i}(eventpoint2:eventpoint2+5.9*60*60*f,:);%6h
    constranit2m=mcherry.correct_NE{i}(eventpoint2:eventpoint2+5.9*60*60*f,:);%6h
    base_DF1=(baseline1g-median(baseline1g))./median(baseline1g)*100;
    base_DF2=(baseline2g-median(baseline2g))./median(baseline2g)*100;
    base_DF3=(baseline1m-median(baseline1m))./median(baseline1m)*100;
    base_DF4=(baseline2m-median(baseline2m))./median(baseline2m)*100;
    constr_DF1=(constranit1g-median(baseline1g))./median(baseline1g)*100;
    constr_DF2=(constranit2g-median(baseline2g))./median(baseline2g)*100;
    constr_DF3=(constranit1m-median(baseline1m))./median(baseline1m)*100;
    constr_DF4=(constranit2m-median(baseline2m))./median(baseline2m)*100;
    gcamp.df_DA{i}=[base_DF1;constr_DF1];
    gcamp.df_NE{i}=[base_DF2;constr_DF2];
    mcherry.df_DA{i}=[base_DF3;constr_DF3];
    mcherry.df_NE{i}=[base_DF4;constr_DF4];
    
    % motion correct
    p = polyfit(mcherry.df_DA{i},gcamp.df_DA{i}, 1);
    G_fit = p(2) + p(1) *  mcherry.df_DA{i};
    F1=gcamp.df_DA{i}- G_fit;
    p = polyfit(mcherry.df_NE{i},gcamp.df_NE{i}, 1);
    G_fit = p(2) + p(1) *  mcherry.df_NE{i};
    F2=gcamp.df_NE{i}- G_fit;
   
    
    gcamp.df_mc_DA{i}=F1;
    gcamp.df_mc_NE{i}=F2;
    total.DF_DA_g(:,i)=gcamp.df_DA{i};
    total.DF_NE_g(:,i)=gcamp.df_NE{i};
    total.DF_DA_m(:,i)=mcherry.df_DA{i};
    total.DF_NE_m(:,i)=mcherry.df_NE{i};
    total.correctDF_DA(:,i)=gcamp.df_mc_DA{i};
    total.correctDF_NE(:,i)=gcamp.df_mc_NE{i};
    name=original(i).name(1:3);
    total.miceID(i,:)=str2num(name);
 
    clearvars -except total gcamp  i mcherry original control event group1 group2 base
end
%% mean/peak df/f plot
%baseline=base min, constr=t min
close all
f=12; t=15;%min
base=2;
%mean
baseline1=mean(total.correctDF_DA((5-base)*60*f+1:5*60*f,:));
baseline2=mean(total.correctDF_NE((5-base)*60*f+1:5*60*f,:));
constr1=mean(total.correctDF_DA(5*60*f+1:(t+5)*60*f,:));
constr2=mean(total.correctDF_NE(5*60*f+1:(t+5)*60*f,:));
meanresult=[baseline1' constr1' baseline2' constr2'];
[h,p1(1)]=ttest(meanresult(:,1),meanresult(:,2));
[h,p1(2)]=ttest(meanresult(:,3),meanresult(:,4));
[h,p1(3)]=ttest(meanresult(:,1),meanresult(:,3));
[h,p1(4)]=ttest(meanresult(:,2),meanresult(:,4));

figure
set(gcf,'unit','centimeters','position',[10 10 8 8 ])
set(gca,'position',[0.2,0.2,0.75,0.75] )
hold on
fillcolor=[255,255,255;255,3,81;255,255,255;2,159 152]/255;
edgecolor=[255,3,81;255,3,81;2,159 152;2,159 152]/255;
position=[1 2 3.5 4.5];
m=mean(meanresult);
err=std(meanresult)./sqrt(size(meanresult,1));
for i=1:4
    bar(position(i),m(i),0.8,'FaceColor', fillcolor(i,:),'EdgeColor',edgecolor(i,:),'LineWidth',1.5)
   [x1(:,i),~]= swarmplot(position(i),meanresult(:,i),50,'k','w');
    errorbar(position(i),m(i),[],err(i),'Color',edgecolor(i,:),'Linewidth',1.5);
end
for j=[1 3]
    for ik=1:length(x1)
        line([x1(ik,j) x1(ik,j+1)],[meanresult(ik,j) meanresult(ik,j+1)],'color','k','linewidth',0.8)
    end
end
xlim([0.5,5]),
xticks([1.5 4]),xticklabels({'DA','NE'});
ylabel(['Mean' '\Delta' 'F/F(%)'])
% ylim([-0.5,1]),yticks([-0.5:0.5:3])
set(gca,'layer','top');
set(gca,'color','none','xcolor',[0 0 0],'ycolor',[0 0 0],'Linewidth',1.6,'FontName','Arial','FontSize',12,'FontWeight','bold','tickdir','out','ticklength',[0.04 0.025]);
y=ylim;
if p1(1)<0.05
    ySig = max(max(meanresult(:,[1 2]))); 
    sigline([1,2],ySig, p1(1));
  
end

if p1(2)<0.05
    ySig = max(max(meanresult(:,[3 4]))); 
    sigline([3,4],ySig, p1(2));
   
end

%peak
baseline1=max(total.correctDF_DA((5-base)*60*f+1:5*60*f,:));
baseline2=max(total.correctDF_NE((5-base)*60*f+1:5*60*f,:));
constr1=max(total.correctDF_DA(5*60*f+1:(t+5)*60*f,:));
constr2=max(total.correctDF_NE(5*60*f+1:(t+5)*60*f,:));
meanresult=[baseline1' constr1' baseline2' constr2'];
figure
set(gcf,'unit','centimeters','position',[10 10 8 8 ])
set(gca,'position',[0.2,0.2,0.75,0.75] )
hold on
fillcolor=[255,255,255;255,3,81;255,255,255;2,159 152]/255;
edgecolor=[255,3,81;255,3,81;2,159 152;2,159 152]/255;
position=[1 2 3.5 4.5];
m=mean(meanresult);
err=std(meanresult)./sqrt(size(meanresult,1));
for i=1:4
    bar(position(i),m(i),0.8,'FaceColor', fillcolor(i,:),'EdgeColor',edgecolor(i,:),'LineWidth',1.5)
    [x1(:,i),~]=swarmplot(position(i),meanresult(:,i),50,'k','w');
    errorbar(position(i),m(i),[],err(i),'Color',edgecolor(i,:),'Linewidth',1.5);
end
for j=[1 3]
    for ik=1:length(x1)
        line([x1(ik,j) x1(ik,j+1)],[meanresult(ik,j) meanresult(ik,j+1)],'color','k','linewidth',0.8)
    end
end
xlim([0.5,5]),
xticks([1.5 4]),xticklabels({'DA','NE'});
ylabel(['Peak' '\Delta' 'F/F(%)'])
% ylim([-0.5,1]),yticks([-0.5:0.5:3])
set(gca,'layer','top');
set(gca,'color','none','xcolor',[0 0 0],'ycolor',[0 0 0],'Linewidth',1.6,'FontName','Arial','FontSize',12,'FontWeight','bold','tickdir','out','ticklength',[0.04 0.025]);

[h,p2(1)]=ttest(meanresult(:,1),meanresult(:,2));
[h,p2(2)]=ttest(meanresult(:,3),meanresult(:,4));
[h,p2(3)]=ttest(meanresult(:,1),meanresult(:,3));
[h,p2(4)]=ttest(meanresult(:,2),meanresult(:,4));
if p2(1)<0.05
    ySig = max(max(meanresult(:,[1 2]))); 
    sigline([1,2],ySig, p2(1));
  
end

if p2(2)<0.05
    ySig = max(max(meanresult(:,[3 4]))); 
    sigline([3,4],ySig, p2(2));
   
end
%% lineplot df/f
close all
f=12; t=15;%min
base=2;
c=[255,3,81;2,159 152]/255;
Lincolor=c;
baseline1=total.correctDF_DA((5-base)*60*f+1:5*60*f,:);
baseline2=total.correctDF_NE((5-base)*60*f+1:5*60*f,:);
constr1=total.correctDF_DA(5*60*f:(t+5)*60*f-1,:);
constr2=total.correctDF_NE(5*60*f:(t+5)*60*f-1,:);
df1=[baseline1;constr1];
df2=[baseline2;constr2];
x=[-1*base:1/(f*60):t];
x(1)=[];
figure
set(gcf,'unit','centimeters','position',[10 10 30 8 ])
set(gca,'position',[0.2,0.2,0.75,0.75] )
hold on
y=df1;
n1=size(y,2);
m1=mean(y,2);
sd1=std(y,1,2);
sem=sd1/sqrt(n1);
shadedErrorBarsub(x,m1,sem,Lincolor(1,:));

y=df2;
n2=size(y,1);
m2=mean(y,2);
sd2=std(y,1,2);
sem=sd2/sqrt(n2);
shadedErrorBarsub(x,m2,sem,Lincolor(2,:));
xlim([-1*base t])
xticks([-2,0:5:15]);
xlabel('Time (min)')
ylabel( ['\Delta' 'F/F(%)'])

ylim([-0.5,1]);yticks([-1:0.5:1]);
line([0 0],[-1,1],'linestyle',':','Linewidth',1.5,'color','k')
set(gca,'color','none','xcolor',[0 0 0],'ycolor',[0 0 0],'Linewidth',1.6,'FontName','Arial','FontSize',12,'FontWeight','bold','tickdir','out')
%% 热图 df/f
close all
f=12; t=15;%min
base=2;
baseline1=total.correctDF_DA((5-base)*60*f+1:5*60*f,:);
baseline2=total.correctDF_NE((5-base)*60*f+1:5*60*f,:);
constr1=total.correctDF_DA(5*60*f:(t+5)*60*f-1,:);
constr2=total.correctDF_NE(5*60*f:(t+5)*60*f-1,:);
df1=[baseline1;constr1];
df2=[baseline2;constr2];
h1=df1';h2=df2';
if length(df1)>3000
    n=20*f;%%downsample
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
    H1=movmean(H1,5,2);
    H2=movmean(H2,5,2);
    H1(:,(base*60*f/n:60*60*12/n:end))=nan;
    H2(:,(base*60*f/n:60*60*12/n:end))=nan;
else
    n=10*f;%%downsample
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
    H1=movmean(H1,5,2);
    H2=movmean(H2,5,2);
    H1(:,(base*60*f/n:60*60*12/n:end))=nan;
    H2(:,(base*60*f/n:60*60*12/n:end))=nan;
end

figure
set(gcf,'unit','centimeters','position',[10 10 20 8 ])
set(gca,'position',[0.1,0.2,0.7,0.7] )
H=heatmap(H1);
H.Colormap=parula;
H.ColorScaling='scaledrows'
H.GridVisible='off'
H.ColorLimits = [0 1];
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
H.ColorScaling='scaledrows';%min=0,max=1
H.GridVisible='off';
% H.ColorbarVisible ='off';
H.ColorLimits = [0 1];
for j=1:length(H.XDisplayLabels)
    H.XDisplayLabels{j}=' ';
end
% for j=1:length(H.YDisplayLabels)
%     H.YDisplayLabels{j}=' ';
% end
H.MissingDataColor = [1 1 1];
H.GridVisible='off';
% H.ColorbarVisible ='off';
%% Z-score plot
f=12; t=15;%min
base=2;
c=[255,3,81;2,159 152]/255;
Lincolor=c;
baseline1=total.correctDF_DA((5-base)*60*f+1:5*60*f,:);
baseline2=total.correctDF_NE((5-base)*60*f+1:5*60*f,:);
constr1=total.correctDF_DA(5*60*f:(t+5)*60*f-1,:);
constr2=total.correctDF_NE(5*60*f:(t+5)*60*f-1,:);
df1=[baseline1;constr1];
df2=[baseline2;constr2];
z1=(df1-mean(baseline1))./std(baseline1);
z2=(df2-mean(baseline2))./std(baseline2);
x=[-1*base:1/(f*60):t];
x(1)=[];
figure
set(gcf,'unit','centimeters','position',[10 10 30 8 ])
set(gca,'position',[0.2,0.2,0.75,0.75] )
hold on
y=z1;
n1=size(y,2);
m1=mean(y,2);
sd1=std(y,1,2);
sem=sd1/sqrt(n1);
shadedErrorBarsub(x,m1,sem,Lincolor(1,:));

y=z2;
n2=size(y,1);
m2=mean(y,2);
sd2=std(y,1,2);
sem=sd2/sqrt(n2);
shadedErrorBarsub(x,m2,sem,Lincolor(2,:));
xlim([-1*base t])
xlabel('Time (min)')
ylabel( 'Z-score')

% ylim([-0.5,1]);yticks([-1:0.5:1]);
line([0 0],[-2,6],'linestyle',':','Linewidth',1.5,'color','k')
set(gca,'color','none','xcolor',[0 0 0],'ycolor',[0 0 0],'Linewidth',1.6,'FontName','Arial','FontSize',12,'FontWeight','bold','tickdir','out')
%% mean/peak Z-score plot plot
%baseline=base min, constr=t min
close all
f=12; t=15;%min
base=2;
%mean
baseline1=total.correctDF_DA((5-base)*60*f+1:5*60*f,:);
baseline2=total.correctDF_NE((5-base)*60*f+1:5*60*f,:);
constr1=total.correctDF_DA(5*60*f+1:(t+5)*60*f,:);
constr2=total.correctDF_NE(5*60*f+1:(t+5)*60*f,:);
df1=[baseline1;constr1];
df2=[baseline2;constr2];
z1=(df1-mean(baseline1))./std(baseline1);
z2=(df2-mean(baseline2))./std(baseline2);
baseline1=mean(z1(1:base*60*f,:));
baseline2=mean(z2(1:base*60*f,:));
constr1=mean(z1(base*60*f+1:end,:));
constr2=mean(z2(base*60*f+1:end,:));


meanresult=[baseline1' constr1' baseline2' constr2'];
[h,p1(1)]=ttest(meanresult(:,1),meanresult(:,2));
[h,p1(2)]=ttest(meanresult(:,3),meanresult(:,4));
[h,p1(3)]=ttest(meanresult(:,1),meanresult(:,3));
[h,p1(4)]=ttest(meanresult(:,2),meanresult(:,4));

figure
set(gcf,'unit','centimeters','position',[10 10 8 8 ])
set(gca,'position',[0.2,0.2,0.75,0.75] )
hold on
fillcolor=[255,255,255;255,3,81;255,255,255;2,159 152]/255;
edgecolor=[255,3,81;255,3,81;2,159 152;2,159 152]/255;
position=[1 2 3.5 4.5];
m=mean(meanresult);
err=std(meanresult)./sqrt(size(meanresult,1));
for i=1:4
    bar(position(i),m(i),0.8,'FaceColor', fillcolor(i,:),'EdgeColor',edgecolor(i,:),'LineWidth',1.5)
   [x1(:,i),~]= swarmplot(position(i),meanresult(:,i),50,'k','w');
    errorbar(position(i),m(i),[],err(i),'Color',edgecolor(i,:),'Linewidth',1.5);
end
for j=[1 3]
    for ik=1:length(x1)
        line([x1(ik,j) x1(ik,j+1)],[meanresult(ik,j) meanresult(ik,j+1)],'color','k','linewidth',0.8)
    end
end
xlim([0.5,5]),
xticks([1.5 4]),xticklabels({'DA','NE'});
ylabel(['Mean z-score'])
% ylim([-0.5,1]),yticks([-0.5:0.5:3])
set(gca,'layer','top');
set(gca,'color','none','xcolor',[0 0 0],'ycolor',[0 0 0],'Linewidth',1.6,'FontName','Arial','FontSize',12,'FontWeight','bold','tickdir','out','ticklength',[0.04 0.025]);
y=ylim;
if p1(1)<0.05
    ySig = max(max(meanresult(:,[1 2]))); 
    sigline([1,2],ySig, p1(1));
  
end

if p1(2)<0.05
    ySig = max(max(meanresult(:,[3 4]))); 
    sigline([3.5,4.5],ySig, p1(2)); 
end
if p2(3)<0.05
    ySig = max(max(meanresult(:,[1 3]))); 
    sigline([1,3.5],ySig, p2(3));  
end

if p2(4)<0.05
    ySig = max(max(meanresult(:,[2 4]))); 
    sigline([2,4.5],ySig+4, p2(4));  
end

%peak
baseline1=total.correctDF_DA((5-base)*60*f+1:5*60*f,:);
baseline2=total.correctDF_NE((5-base)*60*f+1:5*60*f,:);
constr1=total.correctDF_DA(5*60*f+1:(t+5)*60*f,:);
constr2=total.correctDF_NE(5*60*f+1:(t+5)*60*f,:);
df1=[baseline1;constr1];
df2=[baseline2;constr2];
z1=(df1-mean(baseline1))./std(baseline1);
z2=(df2-mean(baseline2))./std(baseline2);
baseline1=max(z1(1:base*60*f,:));
baseline2=max(z2(1:base*60*f,:));
constr1=max(z1(base*60*f+1:end,:));
constr2=max(z2(base*60*f+1:end,:));
meanresult=[baseline1' constr1' baseline2' constr2'];
figure
set(gcf,'unit','centimeters','position',[10 10 8 8 ])
set(gca,'position',[0.2,0.2,0.75,0.75] )
hold on
fillcolor=[255,255,255;255,3,81;255,255,255;2,159 152]/255;
edgecolor=[255,3,81;255,3,81;2,159 152;2,159 152]/255;
position=[1 2 3.5 4.5];
m=mean(meanresult);
err=std(meanresult)./sqrt(size(meanresult,1));
for i=1:4
    bar(position(i),m(i),0.8,'FaceColor', fillcolor(i,:),'EdgeColor',edgecolor(i,:),'LineWidth',1.5)
    [x1(:,i),~]=swarmplot(position(i),meanresult(:,i),50,'k','w');
    errorbar(position(i),m(i),[],err(i),'Color',edgecolor(i,:),'Linewidth',1.5);
end
for j=[1 3]
    for ik=1:length(x1)
        line([x1(ik,j) x1(ik,j+1)],[meanresult(ik,j) meanresult(ik,j+1)],'color','k','linewidth',0.8)
    end
end
xlim([0.5,5]),
xticks([1.5 4]),xticklabels({'DA','NE'});
ylabel(['Peak z-score'])
% ylim([-0.5,1]),yticks([-0.5:0.5:3])
set(gca,'layer','top');
set(gca,'color','none','xcolor',[0 0 0],'ycolor',[0 0 0],'Linewidth',1.6,'FontName','Arial','FontSize',12,'FontWeight','bold','tickdir','out','ticklength',[0.04 0.025]);

[h,p2(1)]=ttest(meanresult(:,1),meanresult(:,2));
[h,p2(2)]=ttest(meanresult(:,3),meanresult(:,4));
[h,p2(3)]=ttest(meanresult(:,1),meanresult(:,3));
[h,p2(4)]=ttest(meanresult(:,2),meanresult(:,4));
if p2(1)<0.05
    ySig = max(max(meanresult(:,[1 2]))); 
    sigline([1,2],ySig, p2(1));  
end

if p2(2)<0.05
    ySig = max(max(meanresult(:,[3 4]))); 
    sigline([3.5,4.5],ySig, p2(2));  
end
if p2(3)<0.05
    ySig = max(max(meanresult(:,[1 3]))); 
    sigline([1,3.5],ySig, p2(3));  
end

if p2(4)<0.05
    ySig = max(max(meanresult(:,[2 4]))); 
    sigline([2,4.5],ySig+4, p2(4));  
end

%% 热图 zsocre
close all
f=12; t=15;%min
base=2;
baseline1=total.correctDF_DA((5-base)*60*f+1:5*60*f,:);
baseline2=total.correctDF_NE((5-base)*60*f+1:5*60*f,:);
constr1=total.correctDF_DA(5*60*f:(t+5)*60*f-1,:);
constr2=total.correctDF_NE(5*60*f:(t+5)*60*f-1,:);
df1=[baseline1;constr1];
df2=[baseline2;constr2];
z1=(df1-mean(baseline1))./std(baseline1);
z2=(df2-mean(baseline2))./std(baseline2);

h1=z1';h2=z2';
if length(df1)>3000
    n=20*f;%%downsample
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
    H1=movmean(H1,5,2);
    H2=movmean(H2,5,2);
    H1(:,(base*60*f/n:60*60*12/n:end))=nan;
    H2(:,(base*60*f/n:60*60*12/n:end))=nan;
else
    n=10*f;%%downsample
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
    H1=movmean(H1,5,2);
    H2=movmean(H2,5,2);
    H1(:,(base*60*f/n:60*60*12/n:end))=nan;
    H2(:,(base*60*f/n:60*60*12/n:end))=nan;
end

figure
set(gcf,'unit','centimeters','position',[10 10 20 8 ])
set(gca,'position',[0.1,0.2,0.7,0.7] )
H=heatmap(H1);
H.Colormap=parula;
% H.ColorScaling='scaledrows'
H.GridVisible='off'
H.ColorLimits = [0 3];
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
H.ColorLimits = [0 3];
for j=1:length(H.XDisplayLabels)
    H.XDisplayLabels{j}=' ';
end
% for j=1:length(H.YDisplayLabels)
%     H.YDisplayLabels{j}=' ';
% end
H.MissingDataColor = [1 1 1];
H.GridVisible='off';
% H.ColorbarVisible ='off';
