fs=44100; %�����Ĳ���Ƶ��Ϊ44.1khz 
duration=1; %¼����ʱ��
fprintf('���������ʼ¼��1:\n'); 
pause; 
fprintf('¼���С���\n'); 
sd1=audiorecorder(fs,16,1); %fs,16,������
recordblocking(sd1,duration); %3��
get(sd1); %��ʾ¼����Ϣ 
% fprintf('�����С���\n'); 
% play(sd1); %����¼��
% pause(3); %��ͣ�����Է���¼��
% fprintf('¼��1������ϡ�\n'); 
y1=getaudiodata(sd1); 
audiowrite('1b.wav',y1,fs);
disp('���ڼ���ο�ģ��Ĳ���...')
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

disp('���ڼ������ģ��Ĳ���...')
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

disp('���ڽ���ģ��ƥ��...')
dist = zeros(1,2);

for i=1:1
for j=1:2
	dist(i,j) = dtw(test(i).mfcc, ref(j).mfcc);
end
end

disp('���ڼ���ƥ����...')
for i=1:1
	[d,j] = min(dist(i,:));
    if j==1
        fprintf('ʶ����Ϊ��yes\n');
    else 
        fprintf('ʶ����Ϊ��no\n' );
    end       
end