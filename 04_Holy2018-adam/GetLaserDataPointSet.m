

function [setOut]=GetLaserDataPointSet(x_pos,y_pos,samples,phiInput,label,MapName,samplingRange)
%*********************************************************************
%*********************************************************************
%*********************************************************************
%函数功能：以输入点为激光探测点，获取周围环境的距离和角度信息
%输入：1.(x_pos,y_pos)，控制激光测距仪探测点坐标，因为所有地图大小为640×480，
%                       x_pos取值范围为(0,640),y_pos的取值范围为(0,480)
%      2.samples，控制采样点数，扫描一圈的采样点数为 samplingRange/180*samples
%      3.phiInput,控制激光测距仪朝向角，0度时为正方向，正方向为水平向右（x轴正向）
%      4.label，可视化激光测距仪采样过程的图片标记
%      5.MapName,控制仿真地图
%      6.samplingRange，控制激光测距仪采样范围，采样范围为(-samplingRange,samplingRange)
%                       samplingRange取值范围(0,180)
%输出：setOut为2×(samplingRange/180*samples)矩阵
%      第i个数据表示为：setOut(:,i)=|range|
%                                   |theta|
%      range,激光测距仪到障碍物的距离
%      theta，当前数据点的角度，theta=i×γ,γ为激光测距仪的角分辨率
%      其中(range;theta)对应某一障碍点的信息
%作者：Shaofeng Wu （由Boss给的程序封装）
%时间：2017.11.19
%邮箱：shaofeng693@126.com
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
theta0 = (-samplingRange:theta_step:samplingRange);         %与前向夹角
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

setOut=[range0;theta0];
