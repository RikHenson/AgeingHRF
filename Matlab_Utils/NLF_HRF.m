function [fit,x1] = NLF_HRF(S);

% Non-linear fitting of HRF based on zero and first-order expansions of
% amplitude and latency (based on code by Darren Price)
%
% yt is the template
% y is the data
% yt is modified to fit y
% pst is the reference timecourse
% t0 = centre of stretch. Points will be stretched outwards from this point 
% GridShift/GridStretch - just for visualisation
% AbsCor - whether to optimise absolue correlation or positive correlation only
% trange = range within which to calculate the fit (2 element vector) [tmin tmax]
% doplot = 3: plot error space, doplot = 2: plot at each iteration, 1:plot at convergence, 0: no plots
% difflim = convergence criterion; 1e-6 gives good fits in a reasonable time. 
%
% Original code for ERPs by Darren Price
% updated for HRFs by rik.henson@mrc-cbu.cam.ac.uk

try
    yt = S.yt;
    y  = S.y;
catch
    error('Need to pass data in S.y and template in S.yt')
end

try    pst = S.pst;         catch   pst = [0:31]+0.5;   end
try    t0  = S.t0;          catch   t0 = 0;             end
try    trange = S.trange;   catch   S.trange = [0 32];  end

try    GridShift = S.GridShift;     catch     GridShift = linspace(-5,5,20);    end
try    GridStrch = S.GridStrch;     catch     GridStrch = linspace(.1,2,20);    end

try    stepsize = S.stepsize;       catch   stepsize = [1 0.1];                 end
try    difflim = S.difflim;         catch   difflim = 0.1;                      end
try    AbsCor = S.AbsCor;           catch   AbsCor = 0;                         end

try    doplot = S.doplot;   catch   doplot = 0;         end

yt = shiftdim(yt);
y = shiftdim(y);

tmin = min(pst); tmax = max(pst);
tindex = 1:length(y);
tindex = tindex(pst >= trange(1) & pst <= trange(2));

cf = zeros(length(GridShift),length(GridStrch));
if doplot > 2
    for sti = 1:length(GridStrch)
        for shi = 1:length(GridShift)
            Shift = GridShift(shi);
            Stretch = GridStrch(sti);
            
            ti = ((pst-t0)/Stretch) + t0 - Shift;
            yi = shiftdim(interp1(pst,y,ti));
            iy = find(~isnan(yi));
            [~,cy] = intersect(tindex,iy);
            % yi(isnan(yi)) = 0;
            %  Cost
            cf(shi,sti) = corr(y(cy),yi(cy));
        end
    end
    
%     if AbsCor
%         [shf, stf] = find(abs(cf) == max(abs(cf(:))));
%     else
%         [shf, stf] = find(cf == max(cf(:)));
%     end
    
    if doplot > 1
        figure(doplot);
        clf
        %subplot(2,1,1)
        pcolor(GridStrch,GridShift,cf);shading flat
        colorbar
        xlabel('Stretch')
        ylabel('Latency')
        caxis([0 1])
    end
end

%%

% Initial Values
% currval = [
%     GridShift(shf) GridStrch(stf)
%     ];

% Random Initial Values
currval = [
    0 1
    ]; % initial values

ind  = 0;
inda = 0;
cfgrad = 10;

% currcf = cf(shf,stf); % current cost function value
% bestcf = currcf;

currcf = 0; % current cost function value
bestcf = -2;

if doplot > 1 
    figure(doplot)
    %subplot(2,1,1)
    hold on
    p1 = plot(currval(2),currval(1),'wo');
    hold off
end

eInc = exp(1i*linspace(0,2*pi,10))';
Increment = [real(eInc)*stepsize(1),imag(eInc)*stepsize(2)];

Nsteps = 0;

yi = shiftdim(interp1(pst,yt,pst));

while cfgrad > difflim
    
    Nsteps = Nsteps+1;
    
    for inci = 1:size(Increment,1)
        % pos dir
        testval = currval+Increment(inci,:);
        Shift = testval(1);
        Stretch = testval(2);
        ti = ((pst-t0)/Stretch) + t0 - Shift;
        yi = shiftdim(interp1(pst,yt,ti));
        iy = find(~isnan(yi));
        yi(isnan(yi)) = 0; % Just for plotting later (doesn't affect corr)
        %  Cost
        [~,cy] = intersect(tindex,iy);
        if length(cy) < length(tindex)/2
            hc(inci,1) = 0;
        else
            hc(inci,1) = corr(y(cy),yi(cy));
        end
        if AbsCor
            dir = find(abs(hc) == max(abs(hc)),1,'first');
        else
            dir = find(hc == max(hc),1,'first');
        end
    end
    
    MaxCor = max(hc(:));
    if AbsCor    
        updatedcf = max(abs(hc(:)));
    else
        updatedcf = MaxCor;
    end
    
    cfgrad = abs(currcf-updatedcf);
    currcf = updatedcf;
    
    % check if caught in loop
    if currcf > bestcf
        bestcf = currcf;
        currval = currval+Increment(dir,:);
        ind = 1;
        inda = inda + 1; % accelleration index
    else
        ind = ind + 1;
        inda = 1;
    end
    
    if ind > 1
        Increment = Increment * 0.75; % decelerate current direction
    end
    
    if inda > 1
        Increment = Increment .* 1; % accelerate current direction (1 = no accelleration)
    end
    
    if doplot > 1 %doplot == 1 || (doplot == 2 && cfgrad <= difflim)
        amp = detrend(y,0)'/detrend(yi,0)';
        % plot time series
        figure(doplot+1)
        %subplot(2,1,2), 
        clf, hold on
%        p2 = plot(pst(tindex)',x(tindex)'./std(x(tindex)),'b-',pst(tindex)',yi(tindex)'./std(yi(tindex)),'r--');
        p2 = plot(pst(tindex)',y(tindex)' - mean(y),'b-','LineWidth',2);
        p2 = plot(pst(tindex)',amp*(yi(tindex)' - mean(yi)),'r-','LineWidth',2);
        if length(tindex) < length(y)
            p2 = plot(pst',y'-mean(y),'b-','LineWidth',1);
            p2 = plot(pst',amp*(yi'-mean(yi)),'r-','LineWidth',1);
        end
        xlabel('real time');ylabel('Amplitude')
        ylim([-3 3])
        if doplot > 1
            set(p1,'XData',currval(2),'YData',currval(1))
        end
        title(['r=',num2str(bestcf,'%2.7f')])
        drawnow
        fprintf('Shift=%3.2f, Stretch=%3.2f, Amp=%3.2f, Baseline=%3.2f\n',currval(1),currval(2),amp,mean(y)-amp*mean(yi))
        pause(0.01)
        
%        if Nsteps == 1, y0 = yi; end
    end
end

y1 = y(tindex)';
yt1 = yi(tindex)';
amp = detrend(y1,0)/detrend(yt1,0); 
fit.amp_off = mean(y1) - amp*mean(yt1);
fit.amp_scl = amp;
fit.fit = amp*yt1 + fit.amp_off;
fit.lat_off = currval(1);
fit.lat_scl = currval(2);
fit.SSE = sum((y1-fit.fit).^2);
fit.R2 = corr(y1',yt1')^2;
fit.Nsteps = Nsteps;

if doplot > 0
    figure(doplot+1),clf, hold on
    p2 = plot(pst(tindex)',y1(tindex)','b-','LineWidth',2);
    p2 = plot(pst(tindex)',fit.fit,'r-','LineWidth',2);
    if length(tindex) < length(y)
        p2 = plot(pst',y1','b-','LineWidth',1);
        p2 = plot(pst',amp*yi'+fit.baseline,'r-','LineWidth',1);
    end
    xlabel('real time');ylabel('Amplitude')
    title(['r=',num2str(bestcf,'%2.7f')])
    fprintf('Shift=%3.2f, Stretch=%3.2f, Amp=%3.2f, Baseline=%3.2f\n',currval(1),currval(2),amp,fit.baseline)
end

% if doplot == 3
%         % plot time series
%         fig(1,'position',[500 500 800 400])
%         clf
%         p2 = plot(pst(tindex)',x(tindex)'./std(x(tindex)),'b-',pst(tindex)',yi(tindex)'./std(yi(tindex)),'r--');
%         xlabel('real time');ylabel('Amplitude')
%         ylim([-4 4])
% end



