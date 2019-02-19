img=imread('image1.jpg');
image(img);
I=imread('f1.png');
set(handles.pushbutton1,'cdata',I);
handles.output = hObject;