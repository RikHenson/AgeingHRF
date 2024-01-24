function err = can_hrf_fit(vP,y,st)
% fits SPM's canonical HRF based on two gamma functions
% vP = parameters: 1-6 are those in spm_hrf; 7 is amplitude scaling
% y  = FIR data
% st = FIR pst
% usage:
%   [p, err, flag]  = fminsearch(@(vP) can_hrf_fit(vP,y,st))

       dt = 0.1;            % resolution of 0.1s
       t  = st/dt;
       
       hrf = spm_hrf(dt,[vP(1:6) 32]);   
       hrf = hrf / max(hrf);
       Y   = hrf(t) * vP(7);             % vP(7) is height 
                  
       err = sum((Y-y).^2); 
       
       figure(2);clf     
       plot(st,y,'bo'), hold on
       plot(st,Y,'r-')
%        pause(0.05)

       
