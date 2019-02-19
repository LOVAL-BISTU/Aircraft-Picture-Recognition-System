function test_example_CNN
load test_2;
load train3000;
train_x = double(reshape(train_x',48,48,3000))/255;
test_x = double(reshape(test_x',48,48,1000))/255;
train_y = double(train_y');
test_y = double(test_y');

%% ex1 Train a 6c-2s-12c-2s Convolutional neural network 
%will run 1 epoch in about 200 second and get around 11% error. 
%With 100 epochs you'll get around 1.2% error

rand('state',0)

cnn.layers = {
    struct('type', 'i') %
    struct('type', 'c', 'outputmaps', 6, 'kernelsize', 5) %convolution layer
    struct('type', 's', 'scale', 2) %sub sampling layer
    struct('type', 'c', 'outputmaps', 12, 'kernelsize',7) %convolution layer
    struct('type', 's', 'scale', 2) %subsampling layer
};


opts.alpha = 0.3;%学习率
opts.batchsize = 20;%批处理数量
opts.numepochs = 1000;%迭代次数

cnn = cnnsetup(cnn, train_x, train_y);
cnn = cnntrain(cnn, train_x, train_y, opts);

[er, bad] = cnntest(cnn, test_x, test_y);%错误率

%plot mean squared error
%figure; plot(cnn.rL);
%assert(er<0.12, 'Too big error');
