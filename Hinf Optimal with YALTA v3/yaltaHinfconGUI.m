% clear all
close all
addpath('YALTA v1.0.1\src');
button = questdlg('Please specify the program GUI type!','GUI Type','Use YALTA & HINFCON','Use HINFCON','Help','Use YALTA & HINFCON');
if strcmp(button,'Use YALTA & HINFCON')
    Yalta_GUI;
elseif strcmp(button,'Use HINFCON')
    hinfconOnly;
else
    winopen('UserGuide_GUI.pdf')
end