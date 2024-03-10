% LLA(WGS-84)����ϵתECEF��Earth-Centered, Earth-Fixed��ֱ������ϵ
%
%  Longitude ����      ����λ �ȣ�
%  Latitude  γ��
%  H         ���θ߶ȣ� ��Ե������ ��   ����λ km
%  LLA = [ Longitude ; Latitude ; H ],Ҳ��Ϊ���󣬾������������Ҫת�����������
%  x�᷽��   ��������ɵ���ָ���㾭����         km
%  y�᷽��   ��x���z�ᾭ��������������õ�
%  z�᷽��   �ɵ���ָ��������
%  X = [ x ; y ; z ] ��λ����


% 2007-5-24

function X = Lla2Ecef( LLA )

re = 6378.137e3;                     % �������뾶
e  = 0.0818191908426214957;        % ����ƫ����
rp = re*sqrt(1-e^2);               % ���򼫰뾶

for i = 1:size(LLA,2)
    longitude = LLA(1,i)/180*pi;      % ����λת���� rad
    latitude  = LLA(2,i)/180*pi;

    gamma = re / sqrt( 1-e^2*sin(latitude)^2 );

    x = ( gamma + LLA(3,i) )*cos(latitude)*cos(longitude);    % ��λ km
    y = ( gamma + LLA(3,i) )*cos(latitude)*sin(longitude);
    z = ( ( 1 - e^2 )*gamma + LLA(3,i) )*sin(latitude);

    X(:,i) = [ x ; y ; z ];
end


