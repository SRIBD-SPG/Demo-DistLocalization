% LLA(WGS-84)坐标系转ECEF（Earth-Centered, Earth-Fixed）直角坐标系
%
%  Longitude 经度      （单位 度）
%  Latitude  纬度
%  H         海拔高度（ 相对地球表面 ）   （单位 km
%  LLA = [ Longitude ; Latitude ; H ],也可为矩阵，矩阵的列数是需要转换的坐标个数
%  x轴方向   赤道面内由地心指向零经度线         km
%  y轴方向   由x轴和z轴经过右手螺旋定理得到
%  z轴方向   由地心指向正北极
%  X = [ x ; y ; z ] 单位：米


% 2007-5-24

function X = Lla2Ecef( LLA )

re = 6378.137e3;                     % 地球赤道半径
e  = 0.0818191908426214957;        % 地球偏心率
rp = re*sqrt(1-e^2);               % 地球极半径

for i = 1:size(LLA,2)
    longitude = LLA(1,i)/180*pi;      % 将单位转化成 rad
    latitude  = LLA(2,i)/180*pi;

    gamma = re / sqrt( 1-e^2*sin(latitude)^2 );

    x = ( gamma + LLA(3,i) )*cos(latitude)*cos(longitude);    % 单位 km
    y = ( gamma + LLA(3,i) )*cos(latitude)*sin(longitude);
    z = ( ( 1 - e^2 )*gamma + LLA(3,i) )*sin(latitude);

    X(:,i) = [ x ; y ; z ];
end


