 %data����Ҫ�˲��ľ���Fsʵ�ʲ����ʣ�f�˲��뾶������Ƶ�ʣ��������˲�
function [output]=filter_lowpass(data,Fs,f)
Mix_Signal=data;
Wc=2*f/Fs; %����Butterworth�˲�����
[b,a]=butter(4,Wc);
Signal_Filter=filter(b,a,Mix_Signal);%�˲�
output=Signal_Filter;
end
