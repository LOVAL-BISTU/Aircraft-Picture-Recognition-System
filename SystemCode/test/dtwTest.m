fs=44100; %声音的采样频率为44.1khz 
duration=1; %录音的时间
fprintf('按任意键开始录音1:\n'); 
pause; 
fprintf('录音中……\n'); 
sd1=audiorecorder(fs,16,1); %fs,16,单声道
recordblocking(sd1,duration); %3秒
get(sd1); %显示录音信息 
% fprintf('放音中……\n'); 
% play(sd1); %播放录音
% pause(3); %暂停程序以放完录音
% fprintf('录音1播放完毕。\n'); 
y1=getaudiodata(sd1); 
audiowrite('1b.wav',y1,fs);
disp('正在计算参考模板的参数...')
for i=1:2
	fname = sprintf('%d.wav',i);
	x=fname;
    [x,fs]=wavread(x);
	[x1 x2] = vad(x);
	m = mfcc(x);
	m = m(x1-2:x2-2,:);
	ref(i).mfcc = m;
   % soundview(x);
end

disp('正在计算测试模板的参数...')
for i=1:1
    %[x,fs]=wavread('E:\\3.wav')
	fname = sprintf('%db.wav',i);
	x=fname;
    [x,fs]=wavread(x);    
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
    else 
        fprintf('识别结果为：no\n' );
    end       
end