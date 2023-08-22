function [xBF] = spm_get_bf(xBF)
% Fill in basis function structure
% FORMAT [xBF] = spm_get_bf(xBF)
%
% xBF.dt      - time bin length {seconds}
% xBF.name    - description of basis functions specified
%               'hrf'
%               'hrf (with time derivative)'
%               'hrf (with time and dispersion derivatives)'
%               '3Gamma'         (initial dip - peak - undershoot; RH)
%               'Fourier set'
%               'Fourier set (Hanning)'
%               'Gamma functions'
%               'Finite Impulse Response'
%               (any other specification will default to 'hrf')
% xBF.length  - window length (seconds)
% xBF.order   - order
% xBF.T       - microtime resolution (for 'hrf*')
%
% xBF.bf      - array of basis functions
%__________________________________________________________________________
%
% spm_get_bf prompts for basis functions to model event or epoch-related
% responses.  The basis functions returned are unitary and orthonormal
% when defined as a function of peri-stimulus time in time-bins.
% It is at this point that the distinction between event and epoch-related 
% responses enters.
%__________________________________________________________________________
% Copyright (C) 1999-2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_get_bf.m 4473 2011-09-08 18:07:45Z guillaume $
 
 
%-Length of time bin
%--------------------------------------------------------------------------
if ~nargin
    str    = 'time bin for basis functions {secs}';
    xBF.dt = spm_input(str,'+1','r',1/16,1);
end
dt = xBF.dt;
 
 
%-Assemble basis functions
%==========================================================================
 
%-Model event-related responses
%--------------------------------------------------------------------------
if ~isfield(xBF,'name')
    spm_input('Hemodynamic Basis functions...',1,'d')
    Ctype = {
        'hrf',...
        'hrf (with time derivative)',...
        'hrf (with time and dispersion derivatives)',...
        '3gamma',...
        '2Gamma (with time derivative)',...
        'Fourier set',...
        'Fourier set (Hanning)',...
        'Gamma functions',...
        'Finite Impulse Response',...
        'User'};
    str   = 'Select basis set';
    Sel   = spm_input(str,2,'m',Ctype);
    xBF.name = Ctype{Sel};
end
 
%-Get order and length parameters
%--------------------------------------------------------------------------
switch xBF.name
 
    case {'3Gamma','Fourier set','Fourier set (Hanning)',...
          'Gamma functions','Finite Impulse Response'}
    %----------------------------------------------------------------------
    try,   l          = xBF.length;
    catch, l          = spm_input('window length {secs}',3,'e',32);
           xBF.length = l;
    end
    try,   h          = xBF.order;
    catch, h          = spm_input('order',4,'e',4);
           xBF.order  = h;
    end
end


%-Create basis functions
%==========================================================================
switch xBF.name
 
    case {'3Gamma'} % Added by RH
        %----------------------------------------------------------------------
        pst   = [0:dt:32]';
        
        P{1} = [6  1];     % Peak
        P{2} = [16 1];     % Post-peak undershoot
        P{3} = [1  0.2];   % Initial undershoot (worthwhile?) (parameters guessed)
                     
        bf = [];
        for b=1:length(P)
            h = spm_Gpdf(pst/dt,P{b}(1)/P{b}(2),dt/P{b}(2));
            %    bf(:,b) = h'/sum(h);
            bf(:,b) = h'/max(h);  % normalise by height rather than area
        end
        
        xBF.length   = size(bf,1)*dt;
        xBF.order    = size(bf,2);
        

    case {'2Gamma (with time derivative)'} % Added by RH
        %----------------------------------------------------------------------
        pst   = [0:dt:32]';
        
        P{1} = [6  1];     % Peak
        P{2} = [16 1];     % Post-peak undershoot
                     
        bf = [];
        for b=1:length(P)
            h = spm_Gpdf(pst/dt,P{b}(1)/P{b}(2),dt/P{b}(2));
            %    bf(:,b) = h'/sum(h);
            bf(:,b) = h'/max(h);  % normalise by height rather than area
        end
        bf(:,2) = -bf(:,2); % make clear undershoot (doesn't matter for fitting)
        
        for b=1:length(P)
            h = spm_Gpdf((pst+1)/dt,P{b}(1)/P{b}(2),dt/P{b}(2));
            h = h/max(h);  % normalise by height rather than area
            D        = bf(:,b) - h;
            bf       = [bf D(:)];
        end
        
        xBF.length   = size(bf,1)*dt;
        xBF.order    = size(bf,2);


    case {'Fourier set','Fourier set (Hanning)'}
    %----------------------------------------------------------------------
    pst   = [0:dt:l]';
    pst   = pst/max(pst);
 
    %-Hanning window
    %----------------------------------------------------------------------
    if strcmp(xBF.name,'Fourier set (Hanning)')
        g = (1 - cos(2*pi*pst))/2;
    else
        g = ones(size(pst));
    end
 
    %-Zeroth and higher Fourier terms
    %----------------------------------------------------------------------
    bf    = g;
    for i = 1:h
        bf = [bf g.*sin(i*2*pi*pst)];
        bf = [bf g.*cos(i*2*pi*pst)];   
    end
 
    case {'Gamma functions'}
    %----------------------------------------------------------------------
    pst   = [0:dt:l]';
    bf    = spm_gamma_bf(pst,h);
 
    case {'Finite Impulse Response'}
    %----------------------------------------------------------------------
    bin   = l/h;
    bf    = kron(eye(h),ones(round(bin/dt),1));
 
    case {'User'} % Hack by RH
    %----------------------------------------------------------------------
    bf    = xBF.bf;
%    warning('Assuming User-specified HRF has values every %d secs',xBF.dt)

    case {'NONE'}
    %----------------------------------------------------------------------
    bf    = 1;
    
otherwise
    %-Microtime resolution
    %----------------------------------------------------------------------
    try
        fMRI_T   = xBF.T;
    catch
        fMRI_T   = spm_get_defaults('stats.fmri.t');
    end
    
    %-Canonical hemodynamic response function
    %----------------------------------------------------------------------
    [bf, p]      = spm_hrf(dt,[],fMRI_T);
 
    %-Add time derivative
    %----------------------------------------------------------------------
    if strfind(xBF.name,'time')

        dp       = 1;
        p(6)     = p(6) + dp;
        D        = (bf(:,1) - spm_hrf(dt,p,fMRI_T))/dp;
        bf       = [bf D(:)];
        p(6)     = p(6) - dp;
 
        %-Add dispersion derivative
        %------------------------------------------------------------------
        if strfind(xBF.name,'dispersion')
 
            dp   = 0.01;
            p(3) = p(3) + dp;
            D    = (bf(:,1) - spm_hrf(dt,p,fMRI_T))/dp;
            bf   = [bf D(:)];
        end
    end
 
    %-Length and order
    %----------------------------------------------------------------------
    xBF.length   = size(bf,1)*dt;
    xBF.order    = size(bf,2);
 
end
 
 
%-Orthogonalise and fill in basis function structure
%--------------------------------------------------------------------------
%xBF.bf = bf;
xBF.bf = spm_orth(bf);
 
 
%==========================================================================
%- S U B - F U N C T I O N S
%==========================================================================

function bf = spm_gamma_bf(u,h)
% Return basis functions (Gamma functions) used for Volterra expansion
% FORMAT bf = spm_gamma_bf(u,h);
% u   - times {seconds}
% h   - order
% bf  - basis functions (mixture of Gammas)
%__________________________________________________________________________
u     = u(:);
bf    = [];
for i = 2:(1 + h)
        m   = 2^i;
        s   = sqrt(m);
        bf  = [bf spm_Gpdf(u,(m/s)^2,m/s^2)];
end


