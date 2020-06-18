

function [setOut]=GetLaserDataPointSet(x_pos,y_pos,samples,phiInput,label,MapName,samplingRange)
%*********************************************************************
%*********************************************************************
%*********************************************************************
%�������ܣ��������Ϊ����̽��㣬��ȡ��Χ�����ľ���ͽǶ���Ϣ
%���룺1.(x_pos,y_pos)�����Ƽ�������̽������꣬��Ϊ���е�ͼ��СΪ640��480��
%                       x_posȡֵ��ΧΪ(0,640),y_pos��ȡֵ��ΧΪ(0,480)
%      2.samples�����Ʋ���������ɨ��һȦ�Ĳ�������Ϊ samplingRange/180*samples
%      3.phiInput,���Ƽ������ǳ���ǣ�0��ʱΪ������������Ϊˮƽ���ң�x������
%      4.label�����ӻ��������ǲ������̵�ͼƬ���
%      5.MapName,���Ʒ����ͼ
%      6.samplingRange�����Ƽ������ǲ�����Χ��������ΧΪ(-samplingRange,samplingRange)
%                       samplingRangeȡֵ��Χ(0,180)
%�����setOutΪ2��(samplingRange/180*samples)����
%      ��i�����ݱ�ʾΪ��setOut(:,i)=|range|
%                                   |theta|
%      range,�������ǵ��ϰ���ľ���
%      theta����ǰ���ݵ�ĽǶȣ�theta=i����,��Ϊ�������ǵĽǷֱ���
%      ����(range;theta)��Ӧĳһ�ϰ������Ϣ
%���ߣ�Shaofeng Wu ����Boss���ĳ����װ��
%ʱ�䣺2017.11.19
%���䣺shaofeng693@126.com
%*********************************************************************
%*********************************************************************
%*********************************************************************



%*********************************************************************
%*********************************************************************
% ��ʼֵ
phi = phiInput; % �������ǳ�����x��нǣ�


%*********************************************************************
%*********************************************************************
% ��ͼ��Ϣ
map = double(imread(MapName));
figure(label); imagesc(map); colormap(gray(256))        
hold on
%set(gca,'xaxislocation','bottom','yaxislocation','left','ydir','reverse') % set origin position 
[yn, xn] = size(map);
max_rang = sqrt(yn*yn+xn*xn);
%*********************************************************************
%*********************************************************************
% ��ʼɨ��

theta_step = 360/samples;
theta0 = (-samplingRange:theta_step:samplingRange);         %��ǰ��н�
range0 = zeros(size(theta0));           %����һ����С��Thetaһ��������
for i = 1:length(theta0)
    alpha = (phi + theta0(i))/180*pi;   %ɨ���߽Ƕ�
    for t = 1:max_rang
        y = y_pos - floor(t*sin(alpha)); %����Ϊ��
        x = x_pos + floor(t*cos(alpha));
        % �ж��Ƿ񳬳���ͼ��Χ
        if y>=1 && y<=yn && x>=1 && x<=xn
            if map(y,x) == 0
                continue;
            end
        end
        % ɨ�赽�߽�
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
