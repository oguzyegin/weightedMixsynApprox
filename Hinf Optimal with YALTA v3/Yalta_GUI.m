function varargout = Yalta_GUI(varargin)
% YALTA_GUI MATLAB code for Yalta_GUI.fig
%      YALTA_GUI, by itself, creates a new YALTA_GUI or raises the existing
%      singleton*.
%
%      H = YALTA_GUI returns the handle to a new YALTA_GUI or the handle to
%      the existing singleton*.
%
%      YALTA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in YALTA_GUI.M with the given input arguments.
%
%      YALTA_GUI('Property','Value',...) creates a new YALTA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Yalta_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Yalta_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Yalta_GUI

% Last Modified by GUIDE v2.5 16-Feb-2017 20:53:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Yalta_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Yalta_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Yalta_GUI is made visible.
function Yalta_GUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Yalta_GUI (see VARARGIN)
% Choose default command line output for Yalta_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
addpath('D:\University\Matlab\Hinf Optimal with YALTA v3\YALTA v1.0.1\src')
movegui(gcf,'center');
global mainGUI; mainGUI = gcf;
% UIWAIT makes Yalta_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Yalta_GUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% W1 numerator and denominator (W1Num and W1Den are polynomials)
function W1_NUM_Callback(hObject, ~, ~)
global W1Num;
W1Num_str = get(hObject,'String');
W1Num_sym = str2sym(W1Num_str);
W1Num = sym2poly(W1Num_sym);

function W1_NUM_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'Enable', 'Inactive');

function W1_DEN_Callback(hObject, eventdata, handles)
global W1Den;
W1Den_str = get(hObject,'String');
W1Den_sym = str2sym(W1Den_str);
W1Den = sym2poly(W1Den_sym);

function W1_DEN_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'Enable', 'Inactive');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% W2 numerator and denominator (W2Num and W2Den are polynomials)
function W2_NUM_Callback(hObject, eventdata, handles)
global W2Num;
W2Num_str = get(hObject,'String');
W2Num_sym = str2sym(W2Num_str);
W2Num = sym2poly(W2Num_sym);

function W2_NUM_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'Enable', 'Inactive');

function W2_DEN_Callback(hObject, eventdata, handles)
global W2Den;
W2Den_str = get(hObject,'String');
W2Den_sym = str2sym(W2Den_str);
W2Den = sym2poly(W2Den_sym);

function W2_DEN_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'Enable', 'Inactive');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Numerator quasi-polynomial qNum -> [(s, s+1) = s^alphaNum*exp(-hNum*s)+(s^alphaNum+1)]
function q_NUM_Callback(hObject, eventdata, handles)
global save_qNum;
syms s;
qNum_str = get(hObject,'String');
Q_NUM = strsplit(qNum_str,',');
nuOfDelays = length(Q_NUM);
delay_split_qNum = s*ones(1,nuOfDelays);
orders = zeros(1,nuOfDelays);
for k=1:length(Q_NUM)
    delay_split_qNum(k)=str2sym(Q_NUM{k}); % 1st index --> most delayed
    orders(k)=length(sym2poly(delay_split_qNum(k))); % (degree+1) of polynomials
end
highestOrder=max(orders);
save_qNum = zeros(nuOfDelays,highestOrder);
for k=1:nuOfDelays
    currentVector=sym2poly(delay_split_qNum(k));
    highOrder=length(currentVector);
    save_qNum(k,highestOrder-highOrder+1:end)=currentVector;
end

function q_NUM_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'Enable', 'Inactive');

% Numerator quasi-polynomial qDen -> [(s, s+1) = s^alphaDen*exp(-hDen*s)+(s^alphaDen+1)]
function q_DEN_Callback(hObject, eventdata, handles)
global save_qDen;
syms s;
qDen_str = get(hObject,'String');
Q_DEN = strsplit(qDen_str,',');
nuOfDelays = length(Q_DEN);
delay_split_qDen = s*ones(1,nuOfDelays);
orders = zeros(1,nuOfDelays);
for k=1:length(Q_DEN)
    delay_split_qDen(k)=str2sym(Q_DEN{k});
    orders(k)=length(sym2poly(delay_split_qDen(k)));
end
highestOrder=max(orders);
save_qDen = zeros(nuOfDelays,highestOrder);
for k=1:nuOfDelays
    currentVector=sym2poly(delay_split_qDen(k));
    highOrder=length(currentVector);
    save_qDen(k,highestOrder-highOrder+1:end)=currentVector;
end

function q_DEN_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'Enable', 'Inactive');


% --- Executes on button press in initialButton.
function initialButton_Callback(hObject, eventdata, handles)
global save_qNum hNum save_qDen hDen ho;
global Mn Md No;
global alphaNum alphaDen;
global alpha mainGUI;

qNum = save_qNum;
qDen = save_qDen;
set(mainGUI, 'HandleVisibility', 'off');
close all;
set(mainGUI, 'HandleVisibility', 'on');

syms s;
rowNu=size(qNum,1);
index=0;
for k=2:rowNu
    if ~isempty(find(qNum(k,:),1))
        index=index+1;
        delayVectorNum(index) = k-1;
    end
end
if ~exist('delayVectorNum')
    delayVectorNum=[1];
    qNum = [qNum(1,:);zeros(1,size(qNum,2)-1) 1e-30];
    adderNum = 1e-30*exp(-hNum*s);
else
    qNum = [qNum(1,:);qNum(delayVectorNum+1,:)];
    adderNum = 0;
end
rowNu=size(qDen,1);
index=0;
for k=2:rowNu
    if ~isempty(find(qDen(k,:),1))
        index=index+1;
        delayVectorDen(index) = k-1;
    end
end
if ~exist('delayVectorDen')
    delayVectorDen=[1];
    qDen = [qDen(1,:);zeros(1,size(qDen,2)-1) 1e-30];
    adderDen = 1e-30*exp(-hDen*s);
else
    qDen = [qDen(1,:);qDen(delayVectorDen+1,:)];
    adderDen = 0;
end
if length(qNum(1,:)) == 1
    uNum = [];
else
    num = delayFrequencyAnalysisMin(qNum,delayVectorNum,alphaNum,hNum);
    uNum = num.UnstablePoles;
end
reversed = 0;
try
    if isnan(uNum)
        qNum2 = qNum(end:-1:1,:);
        for col = 0:(size(qNum2,2)-1)
           qNum2(:,size(qNum2,2)-col) = (-1)^(col)*qNum2(:,size(qNum2,2)-col);
        end
        num = delayFrequencyAnalysisMin(qNum2,delayVectorNum,alphaNum,hNum);
        uNum = num.UnstablePoles;
        reversed = 1;
    else
    end
catch
end
Mn_Num = poly2sym(real(poly(uNum)),s);
Mn_Den = poly2sym(real(poly(-uNum)),s);
% for k=1:length(uNum)
%     unstablePole = round(uNum(k)*1e+4)*1e-4;
%     if k==1
%         Mn_Num = (s-unstablePole);
%         Mn_Den = (s+unstablePole); 
%     else
%         Mn_Num = Mn_Num*(s-unstablePole);
%         Mn_Den = Mn_Den*(s+unstablePole);
%     end
% end
if reversed==1
    Mn_Num2 = 0; Mn_Den2 = 0;
    for k=1:size(qNum,1)
        Mn_Num2 = Mn_Num2+(poly2sym(qNum(k,:),s)*exp(-hNum*(k-1)*s));
        Mn_Den2 = Mn_Den2+(poly2sym(qNum2(k,:),s)*exp(-hNum*(k-1)*s));
    end
    Mn_Num = Mn_Num* Mn_Num2;
    Mn_Den = Mn_Den* Mn_Den2;
end
if ~exist('Mn_Num','var')
   Mn_Num = 1; 
   Mn_Den = 1; 
end
Mn_Num = Mn_Num*exp(-ho*s);
den=delayFrequencyAnalysisMin(qDen,delayVectorDen,alphaDen,hDen,1,1e-3);
den
uNum = den.UnstablePoles;
if isnan(uNum)
    Md_Num = s/s;
    Md_Den = s/s;
    bull = 1;
    Md = Md_Num/Md_Den;
else
    Md_Num = poly2sym(real(poly(uNum)),s);
    Md_Den = poly2sym(real(poly(-uNum)),s);
end

% for k=1:length(uNum)
%     unstablePole = round(uNum(k)*1e+4)*1e-4;
% 	if k==1
%         Md_Num = (s-unstablePole);
%         Md_Den = (s+unstablePole); 
%     else
%         Md_Num = Md_Num*(s-unstablePole);
%         Md_Den = Md_Den*(s+unstablePole);
% 	end
% end
if ~exist('Md_Num')
    Md_Num = s/s;
    Md_Den = s/s;
    bull = 1;
else
    bull = 0;
end
Mn = Mn_Num/Mn_Den;
Mn_str = evalc('pretty(vpa(Mn,4))');
holder=strfind(Mn_str,'---');
if ~isempty(holder)
    newLine = strfind(Mn_str(holder(end):end),'(');
    if isempty(newLine)
        newLine = strfind(Mn_str(holder(end):end),'s');
    end
    Mn_str=[Mn_str(1:holder(1)-3),Mn_str(holder(1)-2:holder(1)-1),Mn_str(holder(1):holder(end)+newLine(1)-2),Mn_str(holder(end)+newLine(1)-1:end)];
else
    newLine = Mn_str;
end
Md = Md_Num/Md_Den;
if bull == 1
    Md_str = 'Md = 1';
else
    Md_str = evalc('sym2zpk(Md)');
end

alpha = uNum;
M = Mn/Md;
for k=1:(length(delayVectorNum)+1)
    if k==1
        Plant_Num = poly2sym(qNum(1,:),s);
    else
        Plant_Num = Plant_Num+poly2sym(qNum(k,:),s)*exp(-delayVectorNum(k-1)*hNum*s);
    end
end
Plant_Num = Plant_Num - adderNum;

for k=1:(length(delayVectorDen)+1)
    if k==1
        Plant_Den = poly2sym(qDen(1,:),s);
    else
        Plant_Den = Plant_Den+poly2sym(qDen(k,:),s)*exp(-delayVectorDen(k-1)*hDen*s);
    end
end
Plant_Den = Plant_Den - adderDen;
Plant = Plant_Num*exp(-ho*s)/Plant_Den;
No = Plant/M;

set(handles.Mn_STR,'FontName','Courier New')
set(handles.Mn_STR,'FontSize',12)
set(handles.Mn_STR,'String',Mn_str);
set(handles.Md_STR,'FontName','Courier New')
set(handles.Md_STR,'FontSize',12)
if bull==0
    set(handles.Md_STR,'String',Md_str(10:(length(Md_str)-42)));
else
    set(handles.Md_STR,'String',Md_str);
end
set(handles.uipanel2,'Visible','On');

% --- Executes on button press in Continue1.
function Continue1_Callback(hObject, eventdata, handles)
global W1Num;
global W1Den;
global W2Num;
global W2Den;
global alpha;
global gPoints;
global minGamma;
global maxGamma;
global threshold;
global Mn;
global gamma_opt;
global betas;

[gamma_opt,xvalue,yvalue,newFunc,betas] = gammaOpt(W1Num,W1Den,W2Num,W2Den,alpha,gPoints,minGamma,maxGamma,threshold,Mn);
plot(xvalue,yvalue); grid on;
xlabel('\gamma'); ylabel('smin(M)'); title('Optimal case');
set(handles.uipanel5,'Visible','On');
set(handles.g_opt,'String',['Optimal Gamma     =     ',num2str(gamma_opt,'%1.9f\n')]);
set(handles.g_opt,'FontSize',14);
set(handles.g_opt,'FontWeight','bold');
figure(3)
plot(xvalue,yvalue);
xlabel('\gamma');
ylabel('smin(R)');
title('Optimal case');
figure(4)
plot(fliplr(xvalue),newFunc);
xlabel('\gamma');
ylabel('-smin(M)');
title('Negative smin(M) versus \gamma');

function lw1_EDIT_Callback(hObject, eventdata, handles)
global wmin;
wmin = str2double(get(hObject,'String'));

function lw1_EDIT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lw2_EDIT_Callback(hObject, eventdata, handles)
global wmax;
wmax = str2double(get(hObject,'String'));

function lw2_EDIT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nuOfPoints_Callback(hObject, eventdata, handles)
global fPoints;
fPoints=str2double(get(hObject,'String'));

function nuOfPoints_CreateFcn(hObject, eventdata, handles)
global fPoints;
fPoints=1000;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minGamma_Callback(hObject, eventdata, handles)
global minGamma;
minGamma = str2double(get(hObject,'String'));

function minGamma_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxGamma_Callback(hObject, eventdata, handles)
global maxGamma;
maxGamma = str2double(get(hObject,'String'));

function maxGamma_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nuOfGamma_Callback(hObject, eventdata, handles)
global gPoints;
gPoints=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function nuOfGamma_CreateFcn(hObject, eventdata, handles)
global gPoints;
gPoints=1000;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function threshold_Callback(hObject, eventdata, handles)
global threshold;
threshold = str2double(get(hObject,'String'));

function threshold_CreateFcn(hObject, eventdata, handles)
global threshold;
threshold = 1e-6; % default value for threshold
set(hObject,'String','1e-6');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function h_NUM_Callback(hObject, eventdata, handles)
global hNum;
entered = get(hObject,'String');
hNum=str2double(entered);
[hNum]=funHandler(hNum,entered);

function h_NUM_CreateFcn(hObject, eventdata, handles)
global hNum; hNum = 1e-4;
set(hObject,'String',num2str(hNum));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function h_DEN_Callback(hObject, eventdata, handles)
global hDen;
entered = get(hObject,'String');
hDen=str2double(entered);
[hDen]=funHandler(hDen,entered);

function h_DEN_CreateFcn(hObject, eventdata, handles)
global hDen; hDen = 1e-4;
set(hObject,'String',num2str(hDen));
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fPower_NUM_Callback(hObject, eventdata, handles)
global alphaNum;
alphaNum=str2double(get(hObject,'String'));

function fPower_NUM_CreateFcn(hObject, eventdata, handles)
global alphaNum;
alphaNum=1;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fPower_DEN_Callback(hObject, eventdata, handles)
global alphaDen;
alphaDen=str2double(get(hObject,'String'));

function fPower_DEN_CreateFcn(hObject, eventdata, handles)
global alphaDen;
alphaDen=1;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function h_PLANT_Callback(hObject, eventdata, handles)
global ho;
entered = get(hObject,'String');
ho=-str2double(entered);
[ho]=funHandler(ho,entered);

function h_PLANT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Works as placeholder
function W2_NUM_ButtonDownFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'On');
w2_num=get(hObject,'String');
if strcmp(w2_num,'W2 Numerator')
   set(hObject,'String',''); 
end
uicontrol(hObject);
function q_NUM_ButtonDownFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'On');
q_num=get(hObject,'String');
if strcmp(q_num,'Plant Numerator')
   set(hObject,'String',''); 
end
uicontrol(hObject);
function q_DEN_ButtonDownFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'On');
q_den=get(hObject,'String');
if strcmp(q_den,'Plant Denominator')
   set(hObject,'String',''); 
end
uicontrol(hObject);
function W1_NUM_ButtonDownFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'On');
w1_num=get(hObject,'String');
if strcmp(w1_num,'W1 Numerator')
   set(hObject,'String',''); 
end
uicontrol(hObject);
function W1_DEN_ButtonDownFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'On');
w1_den=get(hObject,'String');
if strcmp(w1_den,'W1 Denominator')
   set(hObject,'String',''); 
end
uicontrol(hObject);
function W2_DEN_ButtonDownFcn(hObject, eventdata, handles)
set(hObject, 'Enable', 'On');
w2_den=get(hObject,'String');
if strcmp(w2_den,'W2 Denominator')
   set(hObject,'String',''); 
end
uicontrol(hObject);


function controllerCalc_Callback(hObject, eventdata, handles)
run('Approximation.m');

function controllerCalculator(handles)
global W1Num W1Den W2Num W2Den;
global alpha fPoints wmin wmax threshold;
global Mn Md No gamma_opt;
global E F G L Ro Hn hat_G_gamma Ca appHd Hd appNo Copt dInf K_opt;
order_inputs = evalin('base','DESIRED_ORDERS');
[warning, E, G, F, L, Ro, Hn, appHd, Hd, hat_G_gamma, ~, Ca, appNo, Copt]=optController(gamma_opt,W1Num,W1Den,W2Num,W2Den,alpha,threshold,wmin,wmax,fPoints,Mn,Md,No,order_inputs);
if strcmp(warning(1),'C')
    str = evalc('sym2zpk(E)');
    set(handles.E_GAMMA,'FontName','Courier New')
    set(handles.E_GAMMA,'FontSize',8)
    set(handles.E_GAMMA,'String',str(10:length(str)-42))
    str = evalc('sym2zpk(G)');
    set(handles.G_GAMMA,'FontName','Courier New')
    set(handles.G_GAMMA,'FontSize',8)
    set(handles.G_GAMMA,'String',str(10:length(str)-42))
    str = evalc('sym2zpk(F)');
    set(handles.F_GAMMA,'FontName','Courier New')
    set(handles.F_GAMMA,'FontSize',8)
    set(handles.F_GAMMA,'String',str(10:length(str)-42))
    str = evalc('sym2zpk(L)');
    set(handles.L_GAMMA,'FontName','Courier New')
    set(handles.L_GAMMA,'FontSize',8)
    set(handles.L_GAMMA,'String',str(10:length(str)-42))
    str = evalc('sym2zpk(hat_G_gamma)');
    set(handles.G_hat_GAMMA,'FontName','Courier New')
    set(handles.G_hat_GAMMA,'FontSize',8)
    set(handles.G_hat_GAMMA,'String',str(10:length(str)-42))
    str = evalc('sym2zpk(Ro)');
    set(handles.Ro_GAMMA,'FontName','Courier New')
    set(handles.Ro_GAMMA,'FontSize',8)
    set(handles.Ro_GAMMA,'String',str(10:length(str)-42))
    str = evalc('sym2zpk(Hn)');
    set(handles.Hn_GAMMA,'FontName','Courier New')
    set(handles.Hn_GAMMA,'FontSize',8)
    if Hn~=0
        set(handles.Hn_GAMMA,'String',str(10:length(str)-42))
    else
        set(handles.Hn_GAMMA,'String',str(8:14));
    end
    Co = 1/(1+Hn);
    str = evalc('sym2zpk(Co)');
    set(handles.Co_GAMMA,'FontName','Courier New')
    set(handles.Co_GAMMA,'FontSize',8)
    if Co~=1
        set(handles.Co_GAMMA,'String',str(10:length(str)-42))
    else
        set(handles.Co_GAMMA,'String',str(8:14));
    end
    str = evalc('sym2zpk(appHd)');
    set(handles.Hd_GAMMA,'FontName','Courier New')
    set(handles.Hd_GAMMA,'FontSize',8)
    set(handles.Hd_GAMMA,'String',str(10:length(str)-42))
    str = evalc('sym2zpk(appNo)');
    set(handles.No_GAMMA,'FontName','Courier New')
    set(handles.No_GAMMA,'FontSize',8)
    set(handles.No_GAMMA,'String',str(10:length(str)-42))
    set(handles.uipanel14,'Visible','on');
    set(handles.uipanel12,'Visible','on');
    app_Copt = sym2tf(Ca);
    nume = real(app_Copt.num{1});
    deno = real(app_Copt.den{1});
    app_Copt = tf(nume,deno);
    Ca = tf2sym(app_Copt);
    app_Copt = sym2zpk(Ca);
    app_Copt.DisplayFormat = 'frequency'
    format long
    s= 1e+45;
    dInf = eval(Ro)/eval(L)/gamma_opt;
    K_opt = 1/L;
else
    errordlg(warning,'Error')
end

function import_all_variables(handles)
global wmin wmax fpoints minGamma maxGamma gpoints threshold;
global save_qNum alphaNum hNum save_qDen alphaDen hDen ho W1Num W1Den W2Num W2Den;
% try
    inputs = evalin('base','INITIAL_VARIABLES');
    save_qNum = take_struct_field(inputs,'q_numerator',save_qNum,handles.q_NUM,2); 
    save_qDen = take_struct_field(inputs,'q_denominator',save_qDen,handles.q_DEN,2); 
    ho = take_struct_field(inputs,'plant_delay',ho,handles.h_PLANT,0); 
    alphaNum = take_struct_field(inputs,'plant_num_fractional_power',alphaNum,handles.fPower_NUM,0); 
    alphaDen = take_struct_field(inputs,'plant_den_fractional_power',alphaDen,handles.fPower_DEN,0); 
    hNum = take_struct_field(inputs,'plant_num_tau_YALTA',hNum,handles.h_NUM,0); 
    hDen = take_struct_field(inputs,'plant_den_tau_YALTA',hDen,handles.h_DEN,0); 
    W1Num = take_struct_field(inputs,'W1_numerator',W1Num,handles.W1_NUM,1); 
    W1Den = take_struct_field(inputs,'W1_denominator',W1Den,handles.W1_DEN,1); 
    W2Num = take_struct_field(inputs,'W2_numerator',W2Num,handles.W2_NUM,1); 
    W2Den = take_struct_field(inputs,'W2_denominator',W2Den,handles.W2_DEN,1); 
    wmin = take_struct_field(inputs,'log_w_min',wmin,handles.lw1_EDIT,0); 
    wmax = take_struct_field(inputs,'log_w_max',wmax,handles.lw2_EDIT,0); 
    fpoints = take_struct_field(inputs,'Bode_plot_sample_no',fpoints,handles.nuOfPoints,0); 
    minGamma = take_struct_field(inputs,'gamma_min',minGamma,handles.minGamma,0); 
    maxGamma = take_struct_field(inputs,'gamma_max',maxGamma,handles.maxGamma,0); 
    gpoints = take_struct_field(inputs,'Gamma_sample_no',gpoints,handles.nuOfGamma,0); 
    threshold = take_struct_field(inputs,'threshold',threshold,handles.threshold,0); 
% catch err
%     err
%     disp('You need to define a structure called INITIAL_VARIABLES with specific fields defined in the user guide.')
% end
    
    
function stru1 = export_all_variables(myBool)
global wmin wmax fpoints minGamma maxGamma gpoints threshold;
global save_qNum alphaNum hNum save_qDen alphaDen hDen ho W1Num W1Den W2Num W2Den Mn Md No;
global alpha gamma_opt betas;
global E F G L Ro Hn hat_G_gamma Ca appHd Hd appNo Copt dInf K_opt;
stru1 = struct;
stru1 = create_struct_field(stru1,'q_numerator',save_qNum); 
stru1 = create_struct_field(stru1,'q_denominator',save_qDen); 
stru1 = create_struct_field(stru1,'plant_delay',ho); 
stru1 = create_struct_field(stru1,'plant_num_fractional_power',alphaNum); 
stru1 = create_struct_field(stru1,'plant_den_fractional_power',alphaDen); 
stru1 = create_struct_field(stru1,'plant_num_tau_YALTA',hNum); 
stru1 = create_struct_field(stru1,'plant_den_tau_YALTA',hDen); 
stru1 = create_struct_field(stru1,'W1_numerator',W1Num); 
stru1 = create_struct_field(stru1,'W1_denominator',W1Den); 
stru1 = create_struct_field(stru1,'W2_numerator',W2Num); 
stru1 = create_struct_field(stru1,'W2_denominator',W2Den); 
stru1 = create_struct_field(stru1,'log_w_min',wmin); 
stru1 = create_struct_field(stru1,'log_w_max',wmax); 
stru1 = create_struct_field(stru1,'Bode_plot_sample_no',fpoints); 
stru1 = create_struct_field(stru1,'gamma_min',minGamma); 
stru1 = create_struct_field(stru1,'gamma_max',maxGamma); 
stru1 = create_struct_field(stru1,'Gamma_sample_no',gpoints); 
stru1 = create_struct_field(stru1,'threshold',threshold); 
if myBool==0
    assignin('base','INITIAL_VARIABLES',stru1);
else
    stru1 = create_struct_field(stru1,'Mn',Mn); 
    stru1 = create_struct_field(stru1,'Md',Md); 
    stru1 = create_struct_field(stru1,'No',No); 
    stru1 = create_struct_field(stru1,'Noa',appNo);
    stru1 = create_struct_field(stru1,'unstable_poles',alpha);
    stru1 = create_struct_field(stru1,'gamma_opt',gamma_opt);
    stru1 = create_struct_field(stru1,'zeros_of_E',betas);
    stru1 = create_struct_field(stru1,'E',E);
    stru1 = create_struct_field(stru1,'G',G);
    stru1 = create_struct_field(stru1,'F',F);
    stru1 = create_struct_field(stru1,'hat_G',hat_G_gamma);
    stru1 = create_struct_field(stru1,'L',L);
    stru1 = create_struct_field(stru1,'Kopt',K_opt);
    stru1 = create_struct_field(stru1,'Ro',Ro);
    stru1 = create_struct_field(stru1,'d_infty',dInf);
    stru1 = create_struct_field(stru1,'Hn',Hn);
    stru1 = create_struct_field(stru1,'Hd',Hd);
    stru1 = create_struct_field(stru1,'Hda',appHd);
    stru1 = create_struct_field(stru1,'Copt',Copt);
    stru1 = create_struct_field(stru1,'Ca',Ca);
    assignin('base','ALL_VARIABLES',stru1);
end

function structure = create_struct_field(structure,fieldName,variable)
if ~isempty(variable)
    structure = setfield(structure, fieldName, variable);
end

function global_var = take_struct_field(structure,field,global_var,gui_obj,type)
syms s;
if isfield(structure, field)
    global_var = getfield(structure,field);
    if(type==0) % double
        if strcmp(field,'plant_delay')
            set(gui_obj,'String',['-',num2str(global_var)]);
        else
            set(gui_obj,'String',num2str(global_var));
        end
    elseif(type==1) % polynomial
        set(gui_obj,'String',num2str(char(vpa(poly2sym(global_var,s),2))));
    elseif(type==2) % matrix
        result = '';
        for k=size(global_var,1):-1:1
            result = strcat(result,char(vpa(poly2sym(global_var(k,:),s),2)));
            if k~=1
                result = strcat(result,', ');
            end
        end
        set(gui_obj,'String',result);
    end
end


% --------------------------------------------------------------------
function program_options_Callback(hObject, eventdata, handles)
% %


% --------------------------------------------------------------------
function hinfcon_only_Callback(hObject, eventdata, handles)
% hObject    handle to hinfcon_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function import_from_workspace_Callback(hObject, eventdata, handles)
import_all_variables(handles);


% --------------------------------------------------------------------
function export_to_workspace_Callback(hObject, eventdata, handles)
% % 


% --------------------------------------------------------------------
function export_initials_Callback(hObject, eventdata, handles)
export_all_variables(0);


% --------------------------------------------------------------------
function export_all_Callback(hObject, eventdata, handles)
export_all_variables(1);


% --------------------------------------------------------------------
function user_guide_ClickedCallback(hObject, eventdata, handles)
winopen('UserGuide_GUI.pdf')


% --------------------------------------------------------------------
function link_to_src_ClickedCallback(hObject, eventdata, handles)
web('http://www.sciencedirect.com/science/article/pii/S2405896316307376/pdf?md5=14453f144077ac3aa385e435cd797479&amp;pid=1-s2.0-S2405896316307376-main.pdf');


% --------------------------------------------------------------------
function reset_program_ClickedCallback(hObject, eventdata, handles)
clear all;
close all;
run('Yalta_GUI.m');


% --------------------------------------------------------------------
function save_all_ClickedCallback(hObject, eventdata, handles)
folder_name = uigetdir;
h = findall(0,'Type','figure');
[fig_titles,file_names] = original_figure_titles();
for k=1:size(h,1)
    if ~isempty(h(k).Number)
        figure(h(k).Number); h1=get(gca,'title'); fig_title=get(h1,'string');
        if isempty(fig_title)
            fig_axes = gca;
            if(strcmp(fig_axes(1).YLabel.String,'Number of unstable poles') && strcmp(fig_axes(1).XLabel.String,'Delay (s)'))
                savefig(h(k),[folder_name,'\',file_names{2}]);
            end
        end
        for n=1:length(fig_titles)
            if strcmp(fig_title,fig_titles{n})
                savefig(h(k),[folder_name,'\',file_names{n}]);
            elseif strcmp(h(k).Name,fig_titles{n})
                savefig(h(k),[folder_name,'\',file_names{n}]);
            end
        end
    else
        if(strcmp(h(k).Name,fig_titles{1}))
            savefig(h(k),[folder_name,'\',file_names{1}]);
        end
    end
end
ALL_VARIABLES = export_all_variables(1);
save([folder_name,'\','ALL_VARIABLES.mat'],'ALL_VARIABLES');


function [fig_titles,file_names] = original_figure_titles()
fig_titles{1} = 'Finds Hinf optimal controller for neutral/retarded systems'; file_names{1} = 'program_window';
fig_titles{2} = 'Stability Window';  file_names{2} = 'stability_window';
fig_titles{3} = 'Root Loci'; file_names{3} = 'root_loci';
fig_titles{4} = 'Optimal case'; file_names{4} = 'gamma_plot';
fig_titles{5} = 'Bode Plots of Controllers'; file_names{5} = 'bode_controller';
fig_titles{6} = 'Nyquist plot with Copt(jw)'; file_names{6}='nyquist_Copt';
fig_titles{7} = 'Nyquist plot with Ca(jw)'; file_names{7}='nyquist_Ca';
fig_titles{8} = 'Bode Plot of Hd(s)'; file_names{8}='bode_approximation_Hd';
fig_titles{9} = 'Bode Plot of No(s)'; file_names{9}='bode_approximation_No';
fig_titles{10} = 'Performance plot'; file_names{10}='performance_plot';

