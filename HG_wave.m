classdef HG_wave
    %Class of HG type wave
    %   Properties include m and n (doubles), intensity and phase (square
    %   matrix Z plots)
    
    properties (SetAccess = private)
        Xdim
        Ydim
        m
        n
        eField
        phase
        phaseC
        intensity
        intensityC
        pplate
    end
    
    methods
        function obj = HG_wave(m,n,X,Y)
            % initialize args1,2 to m,n and calculate phase/intensity
            obj.m = double(m);
            obj.n = double(n);
            obj.Xdim = size(X,2);
            obj.Ydim = size(Y,1);
            obj.eField = hgForm(X,Y,obj.m,obj.n);
            obj.pplate = hgFormpp(X,Y,obj.m,obj.n);
            obj.intensity = abs(obj.eField).^2;
            obj.phase = obj.MidToZero(angle(obj.eField));
            obj = obj.pruneFuncs(obj);
        end
        
        function newHG = resetDims(HG,X,Y)  % In main script HG4=HG2.resetDims(X,Y) 
            %resetDims(o,x,y)
            %   HG = Wave.resetDims(x,y) Recalculates an HG wave with the
            %   newly parameterized x and y, maintaining m and n. 
            %   [PS: Just create a new instance]
            newHG = HG.HG_wave(HG.m,HG.n,X,Y);
        end
        
        function out = Chop(obj,varargin)
            %CHOP(o,var{1,2,3}) 
            %   F = Wave.Chop(M) Replaces all values of M below the cutoff
            %   (max(M)/e^4) with NaN entries.
            %
            %   F = Wave.Chop(M,D,filler) Creates an F of the values in D
            %   where M is not NaN and fills the rest with filler (either
            %   zeros or NaN)
            if nargin == 2      %single function to chop below cutoff
                cutoff = max(varargin{1}(:))/(exp(4));
                varargin{1}(varargin{1}<cutoff) = NaN;
                out = varargin{1};
            elseif nargin == 4  %two functions and a background filler
                if strcmp(varargin{3}, 'zeros')
                    out = zeros(obj.Ydim,obj.Xdim);
                elseif strcmp(varargin{3}, 'NaN')
                    out = NaN(obj.Ydim,obj.Xdim);
                else
                    out = 0;
                    error('Improper "Chop" parameter');
                end
                out(isfinite(varargin{1})) = varargin{2}(isfinite(varargin{1}));
            end
        end
        
        function out = toNaN(~,func)
            %ToNaN(~,matrix)
            %   F = Wave.toNaN(Z) creates an matrix F with same dimensions
            %   1 and 2 of Z where all zeros are replaced by NaN entries.
            out = NaN(size(func,1),size(func,2));
            out(func ~= 0) = func(func ~= 0);
        end
        
        function out = toZeros(~,func)
            %ToZeros(~,matrix)
            %   F = Wave.toZeros(N) creates an matrix F with same
            %   dimensions 1 and 2 of N where all NaN entries are replaced
            %   by zeros.
            out = zeros(size(func,1),size(func,2));
            size(out)
            out(isfinite(func)) = func(isfinite(func));
        end
        
    end
    methods (Static, Hidden)
        function scaled = MidToZero(phase)
            Zmax = max(phase(:));
            Zmin = min(phase(:));
            Zmid = (Zmax-Zmin)/2+Zmin;
            scaled = (phase - Zmid);
        end
        function out = pruneFuncs(obj)
            %Creates intCut and phaseCut (wave.intensityC & ''.phaseC)
            %   versions of intensity and phase.
            obj.intensityC = obj.Chop(obj.intensity);
            obj.phaseC = obj.MidToZero(obj.Chop(obj.intensityC,obj.phase,'NaN'));
            out = obj;
        end
    end
end

