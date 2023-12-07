 %data是想要滤波的矩阵，Fs实际采样率，f滤波半径（截至频率），按列滤波
function [output]=filter_lowpass(data,Fs,f)
Mix_Signal=data;
Wc=2*f/Fs; %生成Butterworth滤波函数
[b,a]=butter(4,Wc);
Signal_Filter=filter(b,a,Mix_Signal);%滤波
output=Signal_Filter;
end
