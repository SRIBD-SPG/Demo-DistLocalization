% ECEF（Earth-Centered, Earth-Fixed）直角坐标系 转 LLA(WGS-84)坐标系
%
%  x轴方向   赤道面内由地心指向零经度线         km
%  y轴方向   由x轴和z轴经过右手螺旋定理得到
%  z轴方向   由地心指向正北极
%  X = [ x ; y ; z ] ,也可为矩阵，矩阵的列数是需要转换的坐标个数
%  Longitude 经度      （单位 度）
%  Latitude  纬度
%  H         海拔高度（ 相对地球表面 ）   (单位 km)
%  LLA = [ Longitude ; Latitude ; H ]

% 2007-5-24
% clc
% close all
% clear all
% X=[-3234.829893;    4513.534327;    3127.077247];%目标
% X=[-3242.598808;    4512.467517;    3120.607467];%干扰源1
% X=[-3257.182066;    4510.450038 ;   3108.394957];%干扰源2
% X=[-3282.076072 ;   4507.177896  ;  3087.026089];%干扰源3
%时刻1
% X=[ -3743.512501;    5355.362864;    3647.402104 ];%卫星1
% X=[-3752.761528;    5321.998182;    3686.523515];%卫星2
%时刻2
% X=[ -3744.422337;    5350.993471;    3652.876892];
% X=[-3753.648298;    5317.593140;    3691.973322 ];
% %时刻3
% X=[-3745.329276;    5346.619153;    3658.348208 ];
% X=[-3754.532167;    5313.183201;    3697.419620];
% %时刻4
% X=[-3746.233317;    5342.239913;    3663.816045];
% X=[-3755.413135;    5308.768368;    3702.862401];
% %时刻5
% X=[-3747.134460;    5337.855755;    3669.280398];
% X=[-3756.291200;    5304.348646;    3708.301663];
% %时刻6
% X=[-3748.032704;    5333.466683;    3674.741263];
% X=[-3757.166365;    5299.924038;    3713.737398];
% %时刻7
% X=[-3748.928050 ;   5329.072700 ;   3680.198634];
% X=[-3758.038627;    5295.494548 ;   3719.169603];
% %时刻8
% X=[ -3749.820498;    5324.673811;    3685.652507];
% X=[-3758.907988 ;   5291.060179 ;   3724.598271];
% %时刻9
% X=[ -3750.710047;    5320.270019;    3691.102874];
% X=[-3759.774447 ;   5286.620936 ;   3730.023399];
% %时刻10
% X=[-3751.596698;    5315.861327 ;   3696.549733];

% X=[-3760.638003;    5282.176822 ;   3735.444979 ];
function LLA =Ecef2Llaa( X )
a = 6378.137e3;                    % 地球赤道半径
e  = 0.0818191908426214957;        % 地球偏心率
b = a*sqrt(1-e^2);                 % 地球极半径

for i = 1:size(X,2)
    p = sqrt( X(1,i)^2 + X(2,i)^2 );
    theta = atan( X(3,i)*a/p/b );
    e_2 = (a^2 - b^2)/b^2;

    latitude = atan( ( X(3,i) + e_2*b*sin(theta)^3 )/( p - e^2*a*cos(theta)^3 ) );

    longitude = atan2( X(2,i),X(1,i) );

    H = p/cos(latitude) - a/sqrt(1-e^2*sin(latitude)^2);
    Longitude = longitude/pi*180;
    Latitude  = latitude/pi*180;

    LLA(:,i) = [ Longitude ; Latitude ; H ];
end


