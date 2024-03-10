function [position, st] = locatDirtPos(S, pos, fs, latlim, lonlim, h, varargin)

%   This is a brief description of the function
%   根据采样信号矩阵，在制定区域网格遍历搜索，返回最大值位置和网格矩阵
%
%   Detailed description of the function goes here, if needed.
%
%   Input Arguments:
%   S -- signal sample matrix, n-by-m matrix, n samples and m nodes
%   pos -- position of each node in wsg coordinate, m-by-3 matrix
%   fs -- sample rate
%   h - height
%   latlim - search range of Latitude
%   lonlim - search range of Longitude
%   ngrid - number of grid
%   varargin - optional input arguments
%
%   Output Arguments:
%   sfft -- averged fft spectrum
%   ...
%
%   Example:
%   [out1, out2] = functionName(in1, in2, 'option1', value1, 'option2', value2);
%
%   See also: relatedFunctionName
%
%   Author: Wenqiang PU
%   Email: wenqiangpu@cuhk.edu.cn



% default settings 

p = inputParser;
addParameter(p,'method', 'GLRT');
addParameter(p,'ngrid', 10);

p.parse(varargin{:});

method = p.Results.method;
ngrid = p.Results.ngrid;


nNode = size(pos,1); % number of node

S = S/norm(S,'fro');

tgtlatv = linspace(latlim(1), latlim(2), ngrid);
tgtlonv = linspace(lonlim(1), lonlim(2), ngrid);

[tgtlatv1, tgtlonv1] = meshgrid(tgtlatv,tgtlonv);

st = zeros(ngrid, ngrid);
delay = zeros(1, nNode);
for m = 1:ngrid
    for n = 1:ngrid
        for i = 1:nNode
            p1 = [tgtlatv(m), tgtlonv(n), h];
            p2 = [pos(i,1), pos(i,2), pos(i,3)];
            delay(i) = distWGS(p1,p2)/physconst("LightSpeed");
        end
        % 根据时延补偿各感知节点信号
        SynS = roundTimeSyn(S, fs, delay, 200);
        switch upper(method)
            case 'GLRT'
                st(m,n) = detcGlrMul(SynS,'r',1);
            case 'DPD'
                st(m,n) = DPDMUL(SynS,'delay',delay,'sampRate',fs,'NStart',1000);
%                   st(m, n) = DPD(S, delay, fs );
        end
    end
end
st_max=max(max(st));
[id_xmax,id_ymax]=find(st==st_max);
tgtlatv1max=tgtlatv1(1,id_xmax).';
tgtlonv1max=tgtlonv1(id_ymax,1);
position = [mean(tgtlatv1max),mean(tgtlonv1max)];
% figure
% imagesc(tgtlatv1(1,:),tgtlonv1(:,1),abs(st.'))
% colorbar
% hold on
% plot(position(1),position(2),'k.','markersize',20)   %标记一个黑色的圆点
% axis([latlim(1) latlim(2) lonlim(1) lonlim(2) ]);
% xlabel('x (m)')
% ylabel('y (m)')
% for mm=1:length(tgtlatv1max)
%     text(tgtlatv1max(mm),tgtlonv1max(mm),['x=',num2str(tgtlatv1max(mm)),newline,'y=',num2str(tgtlonv1max(mm))]);
% end

end