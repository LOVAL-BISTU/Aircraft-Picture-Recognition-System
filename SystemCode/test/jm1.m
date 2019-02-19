function varargout = jm1(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @jm1_OpeningFcn, ...
                   'gui_OutputFcn',  @jm1_OutputFcn, ...
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


function jm1_OpeningFcn(hObject, eventdata, handles, varargin)
h = handles.figure1; %返回其句柄
newIcon = javax.swing.ImageIcon('图标.png')
figFrame = get(h,'JavaFrame'); %取得Figure的JavaFrame。
figFrame.setFigureIcon(newIcon); %修改图标
set(handles.axes2,'visible','on')
img=imread('backimage.jpg');
image(img);
I=imread('选择图片.jpg');
set(handles.pushbutton1,'cdata',I);
II=imread('识别图片.jpg');
set(handles.pushbutton2,'cdata',II);
III=imread('语音反馈.jpg');
set(handles.pushbutton3,'cdata',III);
IV=imread('精确分类.jpg');
set(handles.pushbutton4,'cdata',IV);
handles.output = hObject;
guidata(hObject, handles);
set(handles.axes2,'visible','off')

function varargout = jm1_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

function pushbutton1_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile('*.jpg;*.png;*.bmp', 'Pick an M-file');
testimage=fullfile(pathname,filename);
if isequal(filename,0)
   disp('用户选择取消')
else
   disp(['用户选择：', fullfile(pathname, filename)])
end
A=imread(testimage);
figure,imshow(A);
imdata=imresize(A,[48,48]);
Id=double(imdata);
h_log=fspecial('log',5,0.5);
a=imfilter(Id,h_log,'corr','replicate');
b=uint8(a);
imgray = rgb2gray(b);
testdata=reshape(imgray',1,[]);
imdata=imread('friendimage.jpg');
Id=double(imdata);
h_log=fspecial('log',5,0.5);
a=imfilter(Id,h_log,'corr','replicate');
b=uint8(a);
imgray = rgb2gray(b);
fridata=reshape(imgray',1,[]);
global test_x;
test_x = [testdata;fridata];
test_x = double(reshape(test_x',48,48,2))/255;

function pushbutton2_Callback(hObject, eventdata, handles)
global test_x;
global prenet;
load netshibie;
    net=cnnff(cnn,test_x);
    prenet=net;
    a=net.o;
    a=round(a);
    result=a(1); 
if result ==1
    display('it is a airplane ，right？');
    [y,fs]=audioread('isairplane.wav');
    player=audioplayer(y,fs);
    play(player);
    pause(4);
else
    display('it is not a airplane ，right？');
    [y,fs]=audioread('notairplane.wav');
    player=audioplayer(y,fs);
    play(player);
    pause(4);
end
[y,fs]=audioread('pleasepushbutton.wav');
player=audioplayer(y,fs);
play(player);
pause(4);


function pushbutton3_Callback(hObject, eventdata, handles)
global prenet;
a=round(prenet.o);
result=a(1);
fprintf('录音中……\n'); 
fs=44100;
duration=2;
sd1=audiorecorder(fs,16,1); 
recordblocking(sd1,duration); 
get(sd1);  
y1=getaudiodata(sd1); 
audiowrite('1b.wav',y1,fs);
disp('正在计算参考模板的参数...')
for i=1:2
	fname = sprintf('%d.wav',i);
	x=fname;
    [x,fs]=audioread(x);
	[x1 x2] = vad(x);
	m = mfcc(x);
	m = m(x1-2:x2-2,:);
	ref(i).mfcc = m;
   % soundview(x);
end

disp('正在计算测试模板的参数...')
for i=1:1
	fname = sprintf('%db.wav',i);
	x=fname;
    [x,fs]=audioread(x);    
	[x1 x2] = vad(x);
	m = mfcc(x);
	m = m(x1-2:x2-2,:);
	test(i).mfcc = m;
end

disp('正在进行模板匹配...')
dist = zeros(1,2);

for i=1:1
for j=1:2
	dist(i,j) = dtw(test(i).mfcc, ref(j).mfcc);
end
end

disp('正在计算匹配结果...')
for i=1:1
	[d,j] = min(dist(i,:));
    if j==1
        fprintf('识别结果为：yes\n');
        [y,fs]=audioread('end.wav');
        d=result;
        player=audioplayer(y,fs);
        play(player);
        pause(4);
    else 
        fprintf('识别结果为：no\n' );
        [y,fs]=audioread('retrain.wav');
        
        player=audioplayer(y,fs);
        play(player);
        pause(4);
       d=~result;
    end       
        load netshibie cnn;
        oldffb=cnn.ffb;
        oldffb
        a=round(cnn.o);
        a(1)=d;
        opts.alpha = 0.5;
        opts.batchsize = 1;
        opts.numepochs =500;
        cnn=cnnapplygrads(cnn,opts);
        save netshibie cnn; 
        ffb=cnn.ffb;
        disp(ffb);
        change=ffb-oldffb;
disp('网络更新量为:');
change
        msgbox({'偏置更新量为:',num2str(change)});
end


function pushbutton4_Callback(hObject, eventdata, handles)
global test_x;
 load netfenlei;
    net=cnnff(cnn,test_x);
    b=net.o;
    b=round(b);
    result=b(1);
if result==1;
    display('这是一架螺旋桨飞机');
    [y,fs]=audioread('luoxuanjiang.wav');
    player=audioplayer(y,fs);
    play(player);
    pause(4);
else
    display('这是一架喷气式飞机')
    [y,fs]=audioread('penqishi.wav');
    player=audioplayer(y,fs);
    play(player);
    pause(4);
end

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
delete(gcf);
