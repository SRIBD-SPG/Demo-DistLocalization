function dronePositions = droneFormation(center, numDrones, range, formation, direction)
% 初始化无人机位置
dronePositions = zeros(numDrones, 3);

% 设置无人机的高度
dronePositions(:, 3) = center(3);

% 转换方向为弧度
direction = deg2rad(direction);

switch formation
    case 'Triangle'
        if numDrones ~= 3
            error('Triangle formation requires exactly 3 drones.');
        end
        for i = 1:numDrones
            angle = direction + (i-1) * 2 * pi / 3; % 计算角度
            dronePositions(i, 1) = center(1) + range * cos(angle); % 计算经度
            dronePositions(i, 2) = center(2) + range * sin(angle); % 计算纬度
        end

    case 'Diamond'
        if numDrones ~= 4
            error('Diamond formation requires exactly 4 drones.');
        end
        for i = 1:numDrones
            angle = direction + (i-1) * pi / 2; % 计算角度
            dronePositions(i, 1) = center(1) + range * cos(angle); % 计算经度
            dronePositions(i, 2) = center(2) + range * sin(angle); % 计算纬度
        end

    case 'Arrow'
        if mod(numDrones, 2) == 0
            error('Arrow formation requires an odd number of drones.');
        end
        for i = 1:numDrones
            if i <= (numDrones - 1) / 2
                angle = direction - deg2rad(120); % 箭头左翼
            else
                angle = direction + deg2rad(120); % 箭头右翼
            end
            % r = range * (numDrones + 1 - 2 * abs(i - (numDrones + 1) / 2)) / numDrones; % 计算半径
            r = range *  abs(i - (numDrones + 1) / 2); % 计算半径
            dronePositions(i, 1) = center(1) + r * cos(angle); % 计算经度
            dronePositions(i, 2) = center(2) + r * sin(angle); % 计算纬度
        end

    case 'Line'
        for i = 1:numDrones
            dronePositions(i, 1) = center(1) + ((i-1) / (numDrones-1) - 0.5) * range * cos(direction); % 计算经度
            dronePositions(i, 2) = center(2) + ((i-1) / (numDrones-1) - 0.5) * range * sin(direction); % 计算纬度
        end

    case 'Random'
        for i = 1:numDrones
            angle = 2 * pi * rand; % 生成随机角度
            r = range * sqrt(rand); % 生成随机半径
            dronePositions(i, 1) = center(1) + r * cos(angle); % 计算经度
            dronePositions(i, 2) = center(2) + r * sin(angle); % 计算纬度
        end
        
    case 'Trapezoid'
        if numDrones < 4
            error('Trapezoid formation requires at least 4 drones.');
        end
        
        % 计算前置两个节点
        front_ang = pi/6;
        dronePositions(1, 1) = center(1) + range * cos(direction + front_ang);
        dronePositions(1, 2) = center(2) + range * sin(direction + front_ang);
        dronePositions(2, 1) = center(1) + range * cos(direction - front_ang);
        dronePositions(2, 2) = center(2) + range * sin(direction - front_ang);
        
        % 计算剩余节点
        back_ang = pi/3;
        dronePositions(3:end, 1) = center(1) + range * linspace(cos(direction - pi + back_ang), cos(direction - pi - back_ang), numDrones - 2)';
        dronePositions(3:end, 2) = center(2) + range * linspace(sin(direction - pi + back_ang), sin(direction - pi - back_ang), numDrones - 2)';
   case 'Arc'
        front_ang = pi/18;
        dronePositions(1, 1) = center(1) + range * cos(direction + front_ang);
        dronePositions(1, 2) = center(2) + range * sin(direction + front_ang);
        dronePositions(2, 1) = center(1) + range * cos(direction - front_ang);
        dronePositions(2, 2) = center(2) + range * sin(direction - front_ang);

        % 计算剩余节点
        range1 = 1/2*range;
        back_ang_start = -pi/4;
        back_ang_end = -pi/6*7;
        theta = linspace(back_ang_start, back_ang_end, numDrones - 2);
        temp = sort(range1*rand(1,numDrones-2));
        rho = ones(1,numDrones-2) *range1+temp; 
        [x, y] = pol2cart(theta, rho);
       % 平移离散点，使其以中心点为中心
       dronePositions(3:end, 1) = x.'+center(1);
       dronePositions(3:end, 2) = y.'+center(2);

    otherwise
        error('Unknown formation: %s', formation);

    


end
end