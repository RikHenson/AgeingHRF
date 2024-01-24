function spm_movie(Action,Orient,F,p,nloops,mode)
% Runs a movie of a plane from a series of images
% FORMAT spm_movie
%
% This is a simple command-line utility which you can use if you
% wish.
%
% The movie can be displayed as transverse, coronal or sagittal
% in movie mode or slider mode which allows slide selection of
% the frame
%
% @(#)spm_movie.m  		John Ashburner, Chloe Hutton

if (nargin==0)
        spm_progress_bar('Clear');
	spm_figure('Clear','Interactive');
        Orient= spm_input('Select slice orientation',...
		1,'m','transverse|sagittal|coronal');
	F      = spm_input(Inf,'*.img','select images');
	p      = spm_input('plane #',1);
        smode=spm_input('Use slider to view movie ? (y/n) ',1,'s','n');
        if strcmp(lower(smode(1)),'y')
           nloops=1;
           mode=1;
        else
           nloops = spm_input('# loops',2);
           mode=0;
        end
	set(spm_figure('FindWin','Interactive'),...
		'Name','Movie','Pointer','watch');   
	spm_movie('Load',Orient,F,p,nloops,mode);
	spm_figure('Clear','Interactive');
	return;
end

switch lower(Action), case('load')
%==========================================================
% Open graphics window to display images as loaded

global PRINTSTR
if isempty(PRINTSTR)
   PRINTSTR='print -dpsc2 -append spm.ps';
end;

n = size(F,1);
figwin=spm_figure('FindWin','Graphics');
if isempty(figwin)
   figwin=spm_figure('Create','Graphics','Graphics');
end
spm_figure('Clear','Graphics');

fig = axes('position',[0.25 0.25 0.5 0.5],'Parent',figwin);
axis image;axis off;
disp('Don''t move the Graphics Window off the edge of the');
disp('screen while spm_movie is preparing the sequence.');
disp('This results in a Segmentation Fault.');

% Read first volume to get dimensions
%---------------------------------------------------------
f = deblank(F(1,:));
V = spm_vol(f);

if Orient==2 | Orient==3
% Sagittal or Coronal slice
%==========================================================
% Display progress bar
spm_progress_bar('Init',n,'Preparing','Images Complete');

   if Orient==2
      % Sagittal slice
      %------------------------------------
      C=[0 0 1 0;0 1 0 0;1 0 0 0;0 0 0 1];
      DIM=V.dim([3 2]);
   elseif Orient==3
      % Coronal slice
      %------------------------------------
      C=[0 0 1 0;1 0 0 0;0 -1 0 0;0 0 0 1];
      % bug!!! --> DIM=V.dim([3 2]);
      DIM=V.dim([3 1]);
      p=-p;
   end;

for j=1:n
	set(0,'CurrentFigure',figwin);
	f = deblank(F(j,:));
	V = spm_vol(f);
        C(3,4)=-p;
	img = spm_slice_vol(V,inv(C),DIM,0);
	bar = zeros(size(img,1),2);
	if (j==1)
		s = 64/max(max(img));
 		im=image([bar img*s],'Parent',fig,'Tag','ToDelete');
                set(figwin,'CurrentAxes',fig);
		axis image; axis xy; axis off
		M = moviein(n,fig);
                idata=get(im,'CData');
                idim=size(idata);
                MI=zeros(idim(1)*idim(2),n);
	else
		l = round(size(img,1)*(j-1)/(n-1));
		if l>0, bar(1:l,2)=64; end;
		im=image([bar img*s],'Parent',fig,'Tag','ToDelete');
                set(figwin,'CurrentAxes',fig);
		axis image; axis xy;axis off
	end
	drawnow;
	M(:,j) = getframe(fig);
	idata=get(im,'CData');
        idim=size(idata);
        MI(:,j)=reshape(idata,idim(1)*idim(2),1);
	spm_progress_bar('Set',j);
end

elseif Orient==1
% Transverse slice
%==========================================================
C=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
DIM=V.dim(1:2);
% Display progress bar
spm_progress_bar('Init',n,'Preparing','Images Complete');

for j=1:n
        set(0,'CurrentFigure',figwin);
	f = deblank(F(j,:));
	V = spm_vol(f);
        C(3,4)=-p;
	img = spm_slice_vol(V,inv(C),DIM,0);
	bar = zeros(size(img,2),2);
	if (j==1)
		s = 64/max(max(img));
 		im=image([flipud(bar) rot90(img*s)],...
		         'Parent',fig,'Tag','ToDelete');
                set(figwin,'CurrentAxes',fig);
		axis image; axis off
		M = moviein(n,fig);
                idata=get(im,'CData');
                idim=size(idata);
                MI=zeros(idim(1)*idim(2),n);
	else
                l = round(size(img,2)*(j-1)/(n-1));
		if l>0, bar(1:l,2)=64; end;
		im=image([flipud(bar) rot90(img*s)],...
		'Parent',fig,'Tag','ToDelete');
                set(figwin,'CurrentAxes',fig);
		axis image;axis off
	end
	drawnow;
	M(:,j) = getframe(fig);
        idata=get(im,'CData');
        idim=size(idata);
        MI(:,j)=reshape(idata,idim(1)*idim(2),1);
       	spm_progress_bar('Set',j);
end

end

spm_figure('Clear','Interactive');

% save movie into figure window
if mode
  M = MI;
end
m_struct=struct('movie',M,'filename',F,'dim',idim,'orient',Orient,...
		'nloops',nloops,'mode',mode);
set(figwin, 'UserData', m_struct);
spm_movie('play');

%==========================================================
case('play')
%==========================================================
% play required movie
figwin=spm_figure('FindWin','Graphics');
if isempty(figwin)
   error('The graphics window has died somehow');
end
ms=get(figwin,'Userdata');
[M Orient idim F mode nloops] = deal(ms.movie,ms.orient,ms.dim,...
				     ms.filename,ms.mode,ms.nloops);
n = size(F,1);

% Display movie or slider as required
%----------------------------------------------------------
if mode==0 
  % make axes again
  fig = axes('position',[0.25 0.25 0.5 0.5],'Parent',figwin);
  axis image;axis off;

   h=findobj('Tag','ToDelete');
   if ~isempty(h)
      delete(h);
   end;
   movie(fig,M,nloops)
elseif mode==1
   if n==1 
      min_step=0.5;
      max_step=1;
      min_slide=0;
      max_slide=0.5;
   elseif n==2
      min_step=0.5;
      max_step=1;
      min_slide=0;
      max_slide=n-1;
   elseif n<=10
      min_step=1/(n-1);
      max_step=1;
      min_slide=0;
      max_slide=n-1;
   else
      min_step=1/(n-1);
      max_step=10/(n-1);
      min_slide=0;
      max_slide=n-1;
   end
   h=findobj('Tag','ToDelete');
   if ~isempty(h)
      delete(h);
   end;

   s=uicontrol('Style','slider','Parent',figwin,...
		    'Position',[200 150 200 30],...
		    'Min',min_slide,'Max',max_slide,...
		    'SliderStep',[min_step max_step],...
                    'Callback','spm_movie(''Scroll'');',...
		    'String','Frame number');
   set(0,'CurrentFigure',figwin);
   fig = axes('position',[0.25 0.25 0.5 0.5],'Parent',figwin);
   I=reshape(M(:,1),idim(1),idim(2));
   image(I,'Parent',fig,'Tag','ToDelete');
   axis image; 
   if Orient==2 | Orient==3
      axis xy;
   end;
   axis off;
   frame=sprintf('%s',deblank(F(1,:)));
   t=title(frame,'FontSize',14,'Tag','ToDelete');
   set(t, 'HandleVisibility', 'On');
end

case('scroll')
%==========================================================
global PRINTSTR
if isempty(PRINTSTR)
   PRINTSTR='print -dpsc2 -append spm.ps';
end;

figwin=spm_figure('FindWin','Graphics');
if isempty(figwin)
   figwin=spm_figure('Create','Graphics','Graphics');
end

ms=get(figwin,'Userdata');
[M Orient idim F mode nloops] = deal(ms.movie, ms.orient,ms.dim,...
				     ms.filename,ms.mode,ms.nloops);

g=get(gcbo,'Value');

h=findobj('Tag','ToDelete');
if ~isempty(h)
   delete(h);
end;
fig = axes('position',[0.25 0.25 0.5 0.5],'Parent',figwin);
I=reshape(M(:,floor(g+1)),idim(1),idim(2));
image(I,'Parent',fig,'Tag','ToDelete');
axis image; 
if Orient==2 | Orient==3
   axis xy; 
end;
axis off;
frame=sprintf('%s',deblank(F(floor(g+1),:)));
t=title(frame,'FontSize',14,'Tag','ToDelete');
set(t, 'HandleVisibility', 'On');

%==========================================================
case('save')
%==========================================================
% save required movie
figwin=spm_figure('FindWin','Graphics');
if isempty(figwin)
   error('The graphics window has died somehow');
end
ms=get(figwin,'Userdata');
if nargin < 2
  fname = spm_input('Filename','+1','s');
else
  fname = Orient;
end
[p f e] = fileparts(fname);
M =ms.movie;
save(fullfile(p, [f '.mat']),'M');

%==========================================================
case('saveavi')
%==========================================================
% save required movie as avi file

% this is implemented only in matlab version 6 
if isempty(which('movie2avi'))
  error('No movie2avi function - maybe need matlab 6?')
end
figwin=spm_figure('FindWin','Graphics');
if isempty(figwin)
   error('The graphics window has died somehow');
end
ms=get(figwin,'Userdata');
if nargin < 2
  fname = spm_input('Filename','+1','s');
else
  fname = Orient;
end
[p f e] = fileparts(fname);
movie2avi(ms.movie,fullfile(p, [f '.avi']));

%======================================================================== 
otherwise
%========================================================================
warning('Unknown action string')

%=======================================================================
end;

   
