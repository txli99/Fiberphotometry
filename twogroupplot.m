%%two group plot
clear
cd G:\LTX\fig2-gcamp\CWQ_selectdata\ARS
shock=load('totalDF.mat');
cd G:\LTX\fig2-gcamp\CWQ_selectdata\Ctrl
control=load('totalDF.mat');
y{1}=shock.DF1*100;
y{2}=control.DF1*100;
y{3}=shock.DF2*100;
y{4}=control.DF2*100;
y{5}=shock.DF1*100-shock.DF2*100;%470-410
y{6}=control.DF1*100-control.DF2*100;%470-410
figure('name','470 nm')
set(gcf,'unit','centimeters','position',[10 5 10 8])%7.25等同 graphpad 5cm
set(gca,'position',[0.2,0.2,0.7,0.7] )%设置图片位置，可根据实际微调
hold on
c=[255 160 64;184 86 215];%%在此处更改颜色，为RGB值
Lincolor=c/255;
x=[-20:1/15:60];
n1=size(y{1},2);
m1=mean(y{1},2,'omitnan');
sd1=std(y{1},1,2,'omitnan');
sem=sd1/sqrt(n1);
shadedErrorBarsub(x,m1,sem,Lincolor(1,:));
n1=size(y{2},2);
m1=mean(y{2},2,'omitnan');
sd1=std(y{2},1,2,'omitnan');
sem=sd1/sqrt(n1);
shadedErrorBarsub(x,m1,sem,Lincolor(2,:));
xlim([-20,60]),xticks([-20:20:60]),
xlim([-10,20]),xticks([-10:10:20]),
ylim([-2,7]),yticks([-2:2:6])
line([0,0],[-2,7],'Color',[0,0,0],'Linewidth',1,'linestyle',':')
xlabel('Time (sec)'),ylabel('DF/F (%)')
legend('ARS','Ctrl','Location','northeast');%在这里可以修改标注，如把'Shock'修改为‘ARS'，注意单引号
set(gca,'color','none','xcolor',[0 0 0],'ycolor',[0 0 0],'Linewidth',1.6,'FontName','Arial','FontSize',12,'FontWeight','bold','tickdir','out','ticklength',[0.04 0.025])
figure('name','410 nm')
set(gcf,'unit','centimeters','position',[10 5 10 8])
set(gca,'position',[0.2,0.2,0.7,0.7] )
n1=size(y{3},2);
m1=mean(y{3},2,'omitnan');
sd1=std(y{3},1,2,'omitnan');
sem=sd1/sqrt(n1);
shadedErrorBarsub(x,m1,sem,Lincolor(1,:));
n1=size(y{4},2);
m1=mean(y{4},2,'omitnan');
sd1=std(y{4},1,2,'omitnan');
sem=sd1/sqrt(n1);
shadedErrorBarsub(x,m1,sem,Lincolor(2,:));
xlim([-20,60]),xticks([-20:20:60]),
xlim([-10,20]),xticks([-10:10:20]),
ylim([-2,7]),yticks([-2:2:6])
line([0,0],[-2,7],'Color',[0,0,0],'Linewidth',1,'linestyle',':')
xlabel('Time (sec)'),ylabel('DF/F (%)')
legend('ARS','Ctrl','Location','northeast');%在这里可以修改标注，如把'Shock'修改为‘ARS'，注意单引号
set(gca,'color','none','xcolor',[0 0 0],'ycolor',[0 0 0],'Linewidth',1.6,'FontName','Arial','FontSize',12,'FontWeight','bold','tickdir','out','ticklength',[0.04 0.025])
figure('name','470-410 nm')
set(gcf,'unit','centimeters','position',[10 5 10 8])
set(gca,'position',[0.2,0.2,0.7,0.7] )
n1=size(y{5},2);
m1=mean(y{5},2,'omitnan');
sd1=std(y{5},1,2,'omitnan');
sem=sd1/sqrt(n1);
shadedErrorBarsub(x,m1,sem,Lincolor(1,:));
n1=size(y{6},2);
m1=mean(y{6},2,'omitnan');
sd1=std(y{6},1,2,'omitnan');
sem=sd1/sqrt(n1);
shadedErrorBarsub(x,m1,sem,Lincolor(2,:));
xlim([-20,60]),xticks([-20:20:60]),
xlim([-10,20]),xticks([-10:10:20]),
ylim([-2,7]),yticks([-2:2:6])
line([0,0],[-2,7],'Color',[0,0,0],'Linewidth',1,'linestyle',':')
xlabel('Time (sec)'),ylabel('DF/F (%)')
legend('ARS','Ctrl','Location','northeast');%在这里可以修改标注，如把'Shock'修改为‘ARS'，注意单引号
set(gca,'color','none','xcolor',[0 0 0],'ycolor',[0 0 0],'Linewidth',1.6,'FontName','Arial','FontSize',12,'FontWeight','bold','tickdir','out','ticklength',[0.04 0.025])
%% mean and peak response
y1=y{5};
y2=y{6};
mean1=mean(y1(301:600,:),'omitnan');%0-20s
mean2=mean(y2(301:600,:),'omitnan');%0-20s
peak1=max(y1(301:600,:));
peak2=max(y2(301:600,:));
response_mean=[mean1' [mean2 nan]']
response_peak=[peak1' [peak2 nan]']


