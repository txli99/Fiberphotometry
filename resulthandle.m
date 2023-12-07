%CWQ,20220729
%%total
clear
file1=dir('*original_signals*.csv');
file2=dir('*control_signals*.csv');
data1=[];data2=[];
for i=1:length(file1)
    dataread1=importdata(file1(i).name);
    dataread2=importdata(file2(i).name);
    data1=[data1 dataread1.data(1:4200,1:2)];
    data2=[data2 dataread2.data(1:4200,1:2)];
    
end %导入数据
f=15;%帧率(HZ)
data1=filter_lowpass(data1,f,2);
data2=filter_lowpass(data2,f,2);
baseline_start=110;
baseline_end=120;
session_start=120;
session_end=140;
n1=baseline_end-baseline_start;
n2=session_end-session_start;
baselinetime=[baseline_start*f:baseline_end*f];%基线时间(s)
baseline1=data1(baselinetime,:);%基线数据
baseline2=data2(baselinetime,:);%基线数据
sessiontime=[session_start*f+1:session_end*f];
session1=data1(sessiontime,:);%实验数据
session2=data2(sessiontime,:);%实验数据
Data1=[baseline1;session1];
Data2=[baseline2;session2];

f01=mean(baseline1);
f02=mean(baseline2);
DF1=(Data1-f01)./f01;
DF2=(Data2-f02)./f02;
% a=isnan(DF1(1,:));
% DF1(:,find(a==1))=[];
% DF2(:,find(a==1))=[];
figure
hold on
x=[-n1:1/f:n2];
y1=DF1*100;
m1=mean(y1,2,'omitnan');
sd1=std(y1,1,2,'omitnan');
sem1=sd1/sqrt(size(y1,2));
shadedErrorBar(x,m1,sem1,'lineProps','r')
y2=DF2*100;
m2=mean(y2,2,'omitnan');
sd2=std(y2,1,2,'omitnan');
sem2=sd1/sqrt(size(y2,2));
shadedErrorBar(x,m2,sem2,'lineProps','k')
legend('470nm','410nm','Location','northeast');
ylim([-2,7]);
% line([0,0],[y(1,1),y(1,2)],'color','k')
xlim([-n1,n2]);xticks([-n1:20:n2]);
xlabel('Time(s)'),ylabel('DF/f(%)')
DF1(n1*f-10:n1*f,:)=nan;
figure
a=isnan(DF1(1,:));
DF1(:,find(a==1))=[];
DF2(:,find(a==1))=[];
h=heatmap(DF1'*100);
h.Colormap=parula;
h.ColorLimits=[-5,5];
h.GridVisible='off';
h.MissingDataColor=[1,1,1];
baseline_mean1=mean(DF1(1:n1*f,:),'omitnan');
session_mean1=mean(DF1(n1*f+1:(n1+60)*f,:),'omitnan');
baseline_mean2=mean(DF2(1:n1*f,:),'omitnan');
session_mean2=mean(DF2(n1*f+1:(n1+60)*f,:),'omitnan');
% DF_m=DF1*100-DF2*100;
% baseline_mean=mean(DF_m(1:n1*f,:),'omitnan')'
% session_mean=mean(DF_m(n1*f+1:(n1+60)*f,:),'omitnan')'
%%
zscore1=(DF1-mean(DF1(1:n1*f,:),'omitnan'))./std(DF1(1:n1*f,:),'omitnan');
zscore2=(DF2-mean(DF2(1:n1*f,:),'omitnan'))./std(DF2(1:n1*f,:),'omitnan');
figure
hold on
x=[-n1:1/f:n2];
y1=zscore1;
m1=mean(y1,2,'omitnan');
sd1=std(y1,1,2,'omitnan');
sem1=sd1/sqrt(size(y1,2));
shadedErrorBar(x,m1,sem1,'lineProps','r')
y2=zscore2;
m2=mean(y2,2,'omitnan');
sd2=std(y2,1,2,'omitnan');
sem2=sd1/sqrt(size(y2,2));
shadedErrorBar(x,m2,sem2,'lineProps','k')
legend('470nm','410nm','Location','northeast');
ylim([-2,7]);
baseline_mean1=mean(zscore1(1:n1*f,:),'omitnan');
session_mean1=mean(zscore1(n1*f+1:(n1+10)*f,:),'omitnan');
baseline_mean2=mean(zscore2(1:n1*f,:),'omitnan');
session_mean2=mean(zscore2(n1*f+1:(n1+10)*f,:),'omitnan');
total=[baseline_mean1' baseline_mean2'  session_mean1' session_mean2'];
max_session=max(zscore1(n1*f+1:(n1+10)*f,:))
%% each mice
file1=dir('*original_signals*.csv');
file2=dir('*control_signals*.csv');
data1=[];data2=[];
for i=1:length(file1)
    dataread1=importdata(file1(i).name);
    dataread2=importdata(file2(i).name);
    data1=[data1 dataread1.data(1:4200,1:2)];
    data2=[data2 dataread2.data(1:4200,1:2)];
    
end %导入数据
f=15;%帧率(HZ)
baseline_start=100;
baseline_end=120;
session_start=120;
session_end=180;
n1=baseline_end-baseline_start;
n2=session_end-session_start;
baselinetime=[baseline_start*f:baseline_end*f];%基线时间(s)
baseline1=data1(baselinetime,:);%基线数据
baseline2=data2(baselinetime,:);%基线数据
sessiontime=[session_start*f+1:session_end*f];
session1=data1(sessiontime,:);%实验数据
session2=data2(sessiontime,:);%实验数据
Data1=[baseline1;session1];
Data2=[baseline2;session2];
f01=mean(baseline1);
f02=mean(baseline2);
DF1=(Data1-f01)./f01;
DF2=(Data2-f02)./f02;
figure
hold on
x=[-n1:1/f:n2];
y1=DF1*100;
y2=DF2*100;
k=1;
for i=1:2:size(DF1,2)
    subplot(5,3,k)
    hold on
 oy1=y1(:,i);
 oy2=y1(:,i+1);
 cy1=y2(:,i);
 cy2=y2(:,i+1);
 plot(x,oy1,'r-')
 plot(x,oy2,'m--')
 plot(x,cy1,'k-')
 plot(x,cy2,'c--')
 xlim([-n1,n2])
 ylabel('DF/F (%)'),xlabel('Time(s)')
 legend({'ch1-ori','ch2-ori','ch1-ctr','ch2-ctr'},'Location','northeast')
 subtitle(file1(k).name(6:8))
 k=k+1;
end
