%%%%%%%% Create and present the ellipse using mgl %%%%%%%%%%%%%%%
%
clear all; clc; close all

% Set some parameters we'll use.
screenParams = mglDescribeDisplays;   
params.fixation1 = [0.6 2 1 1 1 0 0]; % (cross) width, linewidth, color, origin
params.fixationDuration = 0.8;
params.ellipseDuration = 0.3;
params.SOA = 0.6;
params.viewDistance = 600;
params.backgroundRGB = [0.5 0.5 0.5];	% RGB of the background.  All values are in the [0,1] range.
params.enableFeedback = 0;
params.feedbackDuration = 1;
params.ellipseRGB = [1 1 1];
params.ellipseLongAxis = 3.5; % the length of ellipse's long axis (degree visual angle)
params.ellipseAspectRatio = 2; % the ratio of long:short axis
params.ellipseAngleDiff = 10; % the angle difference (degree)
params.ellipseAngleReference = 10; 
radiusAxisLong = visangle2stimsize(params.ellipseLongAxis,0,params.viewDistance,...
    screenParams(2).screenSizeMM(1),screenParams(2).screenSizePixel(1));

% Open mgl window
mglOpen;

% clear both buffers to gray
mglClearScreen(0.5);mglFlush;
mglClearScreen(0.5);mglFlush;

% Prevent any key to Matlab terminal
ListenChar(2)

% Clear the keyboard buffer.
mglGetKeyEvent;

% Present the start text and wait for go signal
mglTextSet('Helvetica',50,[1 1 1],0,0,0,0,0,0,0);
mglClearScreen(0.5);  
mglTextDraw('Press any key to start',[0 0]);
mglFlush
mglGetKeyEvent;
keepLooping = true;
while keepLooping
    key = mglGetKeyEvent;
    if ~isempty(key)
        keepLooping = false;
    end
end
mglClearScreen(0.5);  
mglFlush

% Choose the unit of stimulus display
mglVisualAngleCoordinates(params.viewDistance/10,screenParams(2).screenSizeMM/10);

            
% Fixation point
mglClearScreen(0.5); 
mglFixationCross(params.fixation1);
mglFlush;    
startTimeFixation = mglGetSecs;
while (mglGetSecs - startTimeFixation) < params.fixationDuration
end

% First ellipse            
radiusAxisShort1 = radiusAxisLong / params.ellipseAspectRatio;
angleEllipse1 = params.ellipseAngleReference;
ellipseMat1 = ellipseTextureCreate(radiusAxisLong,radiusAxisShort1,0,params.ellipseRGB,params.backgroundRGB);
ellipseTexture1 = mglCreateTexture(ellipseMat1);
mglClearScreen(0.5);  
mglBltTexture(ellipseTexture1,[0 0],0,0,angleEllipse1);   
mglFlush
startTimeEllipse1 = mglGetSecs;
while (mglGetSecs - startTimeEllipse1) < params.ellipseDuration
end
mglClearScreen(0.5);
mglFixationCross(params.fixation1);
mglFlush;
mglDeleteTexture(ellipseTexture1)

% Second ellipse
radiusAxisShort2 = radiusAxisLong / params.ellipseAspectRatio;
angleEllipse2 = angleEllipse1 - params.ellipseAngleDiff;
ellipseMat2 = ellipseTextureCreate(radiusAxisLong,radiusAxisShort2,0,params.ellipseRGB,params.backgroundRGB);
ellipseTexture2 = mglCreateTexture(ellipseMat2);
mglClearScreen(0.5);              
mglBltTexture(ellipseTexture2,[0 0],0,0,angleEllipse2);
while (mglGetSecs - startTimeEllipse1) < params.SOA
end
mglFlush
startTimeellipse2 = mglGetSecs;
while (mglGetSecs - startTimeellipse2) < params.ellipseDuration
end
mglClearScreen(0.5);
mglFixationCross(params.fixation1);
mglFlush;
mglDeleteTexture(ellipseTexture2)

% Press any key to quit
mglGetKeyEvent;
keepLooping = true;
while keepLooping
    key = mglGetKeyEvent;
    if ~isempty(key)
        keepLooping = false;
    end
end
mglClearScreen(0.5);  
mglFlush
ListenChar(0)
mglClose