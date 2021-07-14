classdef LG_wave
    %Class of LG type wave
    %   Properties include l and p (doubles), intensity and phase (square
    %   matrix Z plots)
    
    properties (SetAccess = private)
        Rdim
        PHIdim
        l
        p
        eField
        phase
        phaseC
        intensity
        intensityC
    end
    
    methods
        function obj = LG_wave(l,p,R,PHI)
            %initialize args1,2 to l,p and calculate phase/intensity
            obj.l = double(l);
            obj.p = double(p);
            obj.Rdim = size(R,2);
            obj.PHIdim = size(PHI,1);
            obj.eField = lgForm(R,PHI,obj.l,obj.p);
            obj.intensity = abs(obj.eField).^2;
            obj.phase = 1/pi*obj.MidToZero(unwrap(unwrap(angle(obj.eField),[],2),[],1));
            obj = obj.pruneFuncs(obj);
        end
        
        function newLG = resetDims(LG,R,PHI)
            %resetDims(o,r,phi)
            %   LG = Wave.resetDims(R,PHI) Recalculates an LG wave with the
            %   newly parameterized r and phi, maintaining l and p.
            %   [PS: Just create a new instance]
            newLG = LG.LG_wave(LG.l,Lg.p,R,PHI);
        end
        
        function out = Chop(obj,varargin)
            %CHOP(o,var{1,2,3})
            %   F = Wave.Chop(M) Replaces all values of M below the cutoff
            %   (max(M)/e^5) with NaN entries.
            %
            %   F = Wave.Chop(M,D,filler) Creates an F of the values in D
            %   where M is not NaN and fills the rest with filler (either
            %   zeros or NaN)
            if nargin == 2      %single function to chop below cutoff
                cutoff = max(varargin{1}(:))/(exp(5));
                varargin{1}(varargin{1}<cutoff(1)) = NaN;
                out = varargin{1};
            elseif nargin == 4  %two functions and a background filler
                if strcmp(varargin{3}, 'zeros')
                    out = zeros(obj.PHIdim,obj.Rdim);
                elseif strcmp(varargin{3}, 'NaN')
                    out = NaN(obj.PHIdim,obj.Rdim);
                else
                    %                     error('Improper background filler
                    %                     for chopping %d and
                    %                     %d.',varargin{1});
                    out = 0;
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

