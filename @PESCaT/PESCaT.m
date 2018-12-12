classdef PESCaT < handle
    % PESCAT    CS reconstruction by Projection onto Epigraph Sets and Calibration over Tensors             
    % The PESCAT class holds the data and performs the reconstruction of 
    % randomly undersampled multi-coil culti-acquisition data using a CS 
    % framework by employing a tensor interpolation kernel and projection
    % onto epigraph sets.
    %
    % To inherit an object from PESCAT class, you have to input:
    %   1-undersampled k-space data with size:
    %   PE*PE*slices*number of acquisitions*number of coils
    %   2-sampling mask
    %   3-type of the dataset: 'bSSFP'/'T1'/'ToF'/'phantom'
    %
    %  (c) ICON Lab 2018
    methods
        function obj = PESCaT(varargin)
            if nargin < 2
                ME = MException('PESCaT:insufficientInputs', ...
                    'you must at least give me the subject name and acceleration factor');
                throw(ME)
            else
                obj.parseInputs(varargin{:});
                obj.simParams.isTrained = 0;
            end
        end
        
        function reconPESCaT(obj)
            % RECONPESCAT performs the reconstruction on PESCAT object. 
            % Output is the reconstructed image.
            if ~obj.simParams.isTrained
                obj.train()
            end
            
            obj.optimize();
        end
        function s = saveobj(obj)
            s.simParams.optimParams = obj.optimParams;
            s.simParams.simParams = obj.simParams;
            s.Data.recon = obj.recon;
            s.Data.sampling = obj.sampling;
        end        
    end
    methods(Static)
        function obj = loadobj(s)
            if isstruct(s)
                sampling = s.Data.sampling;
                dataset = s.simParams.simParams.dataset;
                obj = PESCaT(0,sampling,dataset);
                obj.optimParams = s.simParams.optimParams;
                obj.simParams = s.simParams.simParams;
                obj.recon = s.Data.recon;
                obj.sampling = sampling;
            else
                obj = s;
            end
        end
    end
    methods (Hidden = true, Access = private)
        parseInputs(obj,inputParser,varargin);
        train(obj);
        measure(obj,img);
        optimize(obj,varargin);
        outputInfo(obj);
    end
    properties (Access = private)
        tlambda;
        trecon;
        iter;
    end
    properties
        recon;          %reconstruction result
        SURE;           %SURE method parameters, if applicable
        kernel;         %trained tensor kernel used in the reconstruction
        data;           %k-space data input to the reconstruction
        sampling;       %sampling patterns
        optimParams;    %parameters resulting from the optimization 
        simParams;      %parameters given to perform the recon.
    end
end
