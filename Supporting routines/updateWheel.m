%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets response for wheel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%li

function [exitFlag storeFlag] = updateWheel()

global stimulus

exitFlag = 0;
storeFlag = 0;

gain = 0.2;
keyInfo = mglGetKeyEvent();
mouseInfo = mglGetMouseEvent();

switch stimulus.responseType
    case 'wheel'
        x =  mouseInfo.x*mglGetParam('xPixelsToDevice');
        x = x - (mglGetParam('deviceWidth')*0.5);
        x = x * gain;
        stimulus.feedbackRot = mod(stimulus.feedback - x,pi);
    case 'mouse'
        x = mouseInfo.x;%*mglGetParam('xPixelsToDevice');
        y = mouseInfo.y;%*mglGetParam('yPixelsToDevice');
        x = x - 1440/2;%(mglGetParam('screenWidth')*0.5);
        y = y - 900/2;%(mglGetParam('screenHeight')*0.5);
        theta = atan(y/x);
        stimulus.feedbackRot = theta;
end
if ~isempty(keyInfo)
    if strcmp(keyInfo.charCode,'q');
        exitFlag = 1;
    end
    if strcmp(keyInfo.charCode,'z');
        storeFlag = 1;
    end
end