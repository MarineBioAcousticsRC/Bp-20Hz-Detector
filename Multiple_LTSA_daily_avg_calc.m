%To calculate correct daily average when multiple ltsas are used
%[(V1*time)+(V2*(24-time))]/24

format long

% [(3.8235416666666655*(0+(37/60)+(29/3600)))+...
%     (10.397940041555358*(24-(0+(37/60)+(29/3600))))]/24

v1 = 1.4889695550351245;

v2 = 1.6174695863746928; 

time = 18+(16/60)+(15/3600);

[(v1*time)+(v2*(24-time))]/24