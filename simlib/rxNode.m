classdef rxNode < handle
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
% 2023-07-16: Initial Version 0.1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        name            % name of the TxNode
        maxSpeed        % maximum speed of the TxNode in meters per second
        maxAltitude     % maximum altitude of the TxNode in meters
        motionDensity  % motion noise density

        Hist  % history of TxNode, i.e., position, signal, etc
        nCount

        initialSpeed    % initial speed of the TxNode, [Vx, Vy, Vz]
        currentSpeed    % current speed of the TxNode

        % Cartesian coordinates
        initialPosition        % initial location of the TxNode in [x,y,z] format
        currentPosition        % current location of the TxNode in [x,y,z] format



        % Geographic coordinates
        % Change the attribute name from 'currentLocation' to 'currentGeoLocation', 23/09/26 linwj
        initialGeoLocation        % initial location of the RxNode in [latitude, longitude] format
        currentGeoLocation        % current location of the RxNode in [latitude, longitude] format

        Rx  % txsite class

%         currentSigStr % sigal strength in dBm
%         currentSnr % snr in dB
%         currentDelay % delay in second
%         currentFd % Doppler shift in Hz

        sampRate % sampling rate in Hz
    end

    methods
        function obj = rxNode(para)
            % constructor function to create a new TxNode object
            obj.name = para.name;
            if isfield(para, 'maxSpeed')
                obj.maxSpeed = para.maxSpeed;
            else
                obj.maxSpeed = 1e3;
            end
            if isfield(para, 'maxAltitude')
                obj.maxAltitude = para.maxAltitude;
            else
                obj.maxAltitude = 1e4;
            end


            if isfield(para, 'initialGeoLocation')
                obj.initialGeoLocation = para.initialGeoLocation;
            end

            if isfield(para, 'Rx')
                obj.Rx = para.Rx;
            end

            

            obj.initialPosition = para.initialPosition;
            obj.initialSpeed = para.initialSpeed;

            obj.currentPosition = para.initialPosition;
            obj.currentSpeed = para.initialSpeed;

            obj.motionDensity = para.motionDensity;
            

            % History Record
            obj.Hist = struct();
            obj.nCount = 1;
            obj.Hist.Step(obj.nCount).Position = obj.initialPosition;
            obj.Hist.Step(obj.nCount).Speed = obj.initialSpeed;
            obj.Hist.Step(obj.nCount).Time = 0; % time in second
            obj.Hist.Step(obj.nCount).ItervalT = 0; % time interval in second


            % Geographic coordinates
%             location = cartesianToGeoLocation(obj.initialPosition(1), obj.initialPosition(2),obj.initialPosition(3));
%             obj.initialAltitude = location(3);
%             obj.initialGeoLocation = location(1:2);
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
            obj.Hist.Step(obj.nCount ).ItervalT = T;
            obj.Hist.Step(obj.nCount ).Time = obj.Hist.Step(obj.nCount-1).Time + T;
            obj.Hist.Step(obj.nCount ).Position = obj.currentPosition;
            obj.Hist.Step(obj.nCount ).Speed = obj.currentSpeed;

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
        end

        function  setPosition(obj, position) % set current position [x,y,z]
            obj.currentPosition = position(:);
            obj.currentGeoLocation = cartesianToGeoLocation(position(1), position(2), position(3));
        end

        function refreshSite(obj) % update geo location of site class
            obj.Rx.Latitude = obj.currentGeoLocation(1);
            obj.Rx.Longitude = obj.currentGeoLocation(2);
            obj.Rx.AntennaHeight = obj.currentGeoLocation(3);
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
            disp(['Rx Name: ', obj.name]);
            disp(['Current Speed: ', num2str(obj.currentSpeed), ' m/s']);
            disp(['Current Altitude: ', num2str(obj.currentAltitude), ' m']);
            disp(['Current GeoLocation: (', num2str(obj.currentGeoLocation(1)), ', ', num2str(obj.currentGeoLocation(2)), ')']);
        end
    end
end






