classdef txNode < handle
% 
%  Class function of TxNode
%  
%
%   Example:
% 
%
%   NOTE: 
%
%   Author: Wenqiang PU 
%   Email: wenqiangpu@cuhk.edu.cn
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2023-07-09: Initial Version 0.1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        Name            % name of the TxNode
        maxSpeed        % maximum speed of the TxNode in meters per second
        maxAltitude     % maximum altitude of the TxNode in meters
        motionDensity   % motion noise density

        Hist  % history of TxNode, i.e., position, signal, etc
        nCount

        initialSpeed    % initial speed of the TxNode, [Vx, Vy, Vz]
        currentSpeed    % current speed of the TxNode

        % Cartesian coordinates
        initialPosition        % initial location of the TxNode in [x,y,z] format
        currentPosition        % current location of the TxNode in [x,y,z] format
        
        
        % Geographic coordinates 

        initialGeoLocation     % initial location of the TxNode in [latitude, longitude altitude] format
        currentGeoLocation     % current location of the TxNode in [latitude, longitude] format
        
        % basic paramerters of transmit signal
        moduType    % modulation type
        baudRate    % baud rate in Hz
        sampRate    % sample rate in Hz
        
        % transmit pattern
        transPattern    % signal transmiton pattern
        traPatPara       % signal transmition parameters
        
        
        % txsite class 
        Tx  % txsite class
    end

    methods
        function obj = txNode(varargin)
            % constructor function to create a new TxNode object
            p = inputParser;

            addParameter(p, 'Name', 'TxNode'); % name
            addParameter(p, 'maxSpeed', 1e3); % m/s
            addParameter(p, 'maxAltitude', 3e4); % meter 

            addParameter(p, 'initialGeoLocation', [22.806537, 114.250512, 1000]'); % geo location 
            addParameter(p, 'initialPosition', []); %
            addParameter(p, 'initialSpeed', zeros(3,1)); % [vx, vy, vz]

            addParameter(p, 'currentGeoLocation', []); %
            addParameter(p, 'currentPosition', []); %
            addParameter(p, 'currentSpeed', []); %

            addParameter(p, 'motionDensity', 0); % motion noise

            addParameter(p, 'moduType', []); %
            addParameter(p, 'baudRate', []); % 
            addParameter(p, 'sampRate', []); %

            addParameter(p, 'transPattern', []); %
            addParameter(p, 'traPatPara', []); %

            addParameter(p, 'Tx', []); %

            p.parse(varargin{:});

            
            obj.Name = p.Results.Name;
            obj.maxSpeed = p.Results.maxSpeed;
            obj.maxAltitude = p.Results.maxAltitude;
            
            obj.initialGeoLocation = p.Results.initialGeoLocation;
            obj.initialPosition = p.Results.initialPosition;
            obj.initialSpeed = p.Results.initialSpeed;

            obj.currentGeoLocation = p.Results.initialGeoLocation;
            obj.currentPosition = p.Results.initialPosition;
            obj.currentSpeed = p.Results.initialSpeed;

            obj.motionDensity = p.Results.motionDensity;

            obj.moduType = p.Results.moduType;
            obj.baudRate = p.Results.baudRate;
            obj.sampRate = p.Results.sampRate;

            obj.transPattern = p.Results.transPattern;
            obj.traPatPara = p.Results.traPatPara;

            obj.Tx = p.Results.Tx;

            % History Record
            obj.Hist = struct();
            obj.nCount = 1;
            obj.Hist.Step(obj.nCount).Position = obj.initialPosition;
            obj.Hist.Step(obj.nCount).Speed = obj.initialSpeed;
            obj.Hist.Step(obj.nCount).Time = 0; % time in second
            obj.Hist.Step(obj.nCount).ItervalT = 0; % time interval in second


        end

        function moveNode(obj, speed, T)
            % function to move the TxNode to a new location at a given speed and altitude
            if norm(speed) > obj.maxSpeed
                error('Speed exceeds maximum speed for this TxNode');
            end
%             if altitude > obj.maxAltitude
%                 error('Altitude exceeds maximum altitude for this TxNode');
%             end

            obj.currentSpeed = speed;
            obj.currentPosition = obj.currentPosition + T * speed;

            obj.nCount = obj.nCount + 1;
            obj.Hist.Step(obj.nCount).ItervalT = T;
            obj.Hist.Step(obj.nCount).Time = obj.Hist.Step(obj.nCount-1).Time + T;
            obj.Hist.Step(obj.nCount).Position = obj.currentPosition;
            obj.Hist.Step(obj.nCount).Speed = obj.currentSpeed;

%             obj.currentAltitude = altitude;
%             obj.location = location;
        end
        

        function set.currentGeoLocation(obj, location)
            %%% modified by linwj,23/09/26
            % when currentGeoLocation is assigned a new value ,this method will be
            % excuted automaticly 
            obj.currentGeoLocation =location;% assign
            obj.refreshSite(); % update geo location of site class
        end

        function  setGeoLocation(obj, lat, lon, alt) % set current location
            obj.currentGeoLocation = [lat, lon, alt]';
            obj.currentPosition = locationToCartesian(lat, lon, alt);
           
            %%% modified by linwj,23/09/26, 在currentGeoLocation被赋值时，自动执行
            % refreshSite(obj); % update geo location of site class
        end

        function  setPosition(obj, position) % set current position [x,y,z]
            obj.currentPosition = position(:);
            obj.currentGeoLocation = cartesianToGeoLocation(position(1), position(2), position(3));

            %%% modified by linwj,23/09/26, 在currentGeoLocation被赋值时，自动执行
            %refreshSite(obj); % update geo location of site class
        end

        function refreshSite(obj) % update geo location of site class
            obj.Tx.Latitude = obj.currentGeoLocation(1);
            obj.Tx.Longitude = obj.currentGeoLocation(2);
            obj.Tx.AntennaHeight = obj.currentGeoLocation(3);
        end



        function  setSpeed(obj, speed) % set current speed
            obj.currentSpeed = speed(:);
        end



        function  x = getX(obj)
                x = obj.currentPosition(1);
        end

        function  y = getY(obj)
                y = obj.currentPosition(2);
        end

        function  z = getZ(obj)
                z = obj.currentPosition(3);
        end

        function displayNode(obj)
            % function to display the current state of the TxNode
            disp(['Tx Name: ', obj.Name]);
            disp(['Current Speed: ', num2str(obj.currentSpeed), ' m/s']);
            disp(['Current GeoLocation: (', num2str(obj.currentGeoLocation(1)), ', ', num2str(obj.currentGeoLocation(2)), ')']);
        end
    end
end






