

function [setOut]=GetLaserDataPointSet(x_pos,y_pos,samples,phiInput,label,MapName)
%*********************************************************************
%*********************************************************************
%*********************************************************************
%函数功能：以输入点为激光探测点，获取周围环境的距离和角度信息
%输入：激光测距仪探测点坐标 （x_pos,y_pos）
%      扫描一圈的采样点数 samples
%输出：2×1矩阵setOut=|range|
%                     |theta|
%      距离向量 range
%      角度向量 theta
%      其中(range(i),theta(i))对应某一障碍点的信息
%作者：Shaofeng Wu （由Boss给的程序封装）
%时间：2017.11.19
%*********************************************************************
%*********************************************************************
%*********************************************************************



%*********************************************************************
%*********************************************************************
% 初始值
phi = phiInput; % 激光测距仪朝向（与x轴夹角）


%*********************************************************************
%*********************************************************************
% 地图信息
map = double(imread(MapName));
figure(label); imagesc(map); colormap(gray(256))        
hold on
%set(gca,'xaxislocation','bottom','yaxislocation','left','ydir','reverse') % set origin position 
[yn, xn] = size(map);
max_rang = sqrt(yn*yn+xn*xn);
%*********************************************************************
%*********************************************************************
% 开始扫描

theta_step = 360/samples;
theta0 = (-120:theta_step:120);         %与前向夹角
range0 = zeros(size(theta0));           %设置一个大小和Theta一样的向量
for i = 1:length(theta0)
    alpha = (phi + theta0(i))/180*pi;   %扫描线角度
    for t = 1:max_rang
        y = y_pos - floor(t*sin(alpha)); %向下为正
        x = x_pos + floor(t*cos(alpha));
        % 判断是否超出地图范围
        if y>=1 && y<=yn && x>=1 && x<=xn
            if map(y,x) == 0
                continue;
            end
        end
        % 扫描到边界
        range0(i) = t;
        plot([x_pos x],[y_pos y],'-')
        drawnow
        break;
    end
end
hold off


figure(label+1); p = plot(theta0,range0, '-');
set(p,'linewidth',2)

% setOut=[range0;theta0];%输出为角度制，范围为-180~180
setOut=[range0;theta0/180*pi+pi];%输出为弧度制，范围为0~2*pi
