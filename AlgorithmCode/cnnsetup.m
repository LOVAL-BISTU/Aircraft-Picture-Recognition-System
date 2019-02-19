function net = cnnsetup(net, x, y)
    assert(~isOctave() || compare_versions(OCTAVE_VERSION, '3.8.0', '>='), ['Octave 3.8.0 or greater is required for CNNs as there is a bug in convolution in previous versions. See http://savannah.gnu.org/bugs/?39314. Your version is ' myOctaveVersion]);
    inputmaps = 1;
    mapsize = size(squeeze(x(:, :, 1)));%除去size为1的维度然后获取行数列数。
    load dbnmap.mat%修改后添加
    p=1;
    for l = 1 : numel(net.layers)   % layer，每一层。numel函数用于计算数组中满足指定条件的元素个数。若是一幅图像，则numel(A)将给出它的像素数。
        if strcmp(net.layers{l}.type, 's')%c = strcmp(str1,str2)，比较字符串 str1 与 str2 ，若完全相等则返回 1 ，不相等返回 0。
            mapsize = mapsize / net.layers{l}.scale;%这种把采样频率降下来，就是降采样（downsample）。这样做的好处是减少数据样点，也就是减少运算时间，在实时处理时常采用的方法。这里除以scale=2，就是pooling之后图的大小，pooling域之间没有重叠，所以pooling后的图像为25*25.mapsize保存的是上一层的size
            assert(all(floor(mapsize)==mapsize), ['Layer ' num2str(l) ' size must be integer. Actual: ' num2str(mapsize)]);
            for j = 1 : inputmaps
                net.layers{l}.b{j} = 0;
            end
        end
        if strcmp(net.layers{l}.type, 'c')%如果是卷积层。
            mapsize = mapsize - net.layers{l}.kernelsize + 1;%得到featuremap的size，大小为上层大小-卷积核大小+1(即50-7+1=44）（22-11+1）。featruemap是每种滤波器卷积得到不同特征的放映。
            fan_out = net.layers{l}.outputmaps * net.layers{l}.kernelsize ^ 2;%得到该层所有连接的数量=卷积核数*卷积核size=6*（7*7），12（11*11）
            for j = 1 : net.layers{l}.outputmaps  %  output map 卷积核数*卷积核size
                fan_in = inputmaps * net.layers{l}.kernelsize ^ 2;     %本层的一个输出map所对应的所有卷积核，包含的权重神经元的总数=1*7*7，6*11*11
                for i = 1 : inputmaps      %  input map%原CNN为net.layers{l}.k{i}{j} = (rand(net.layers{l}.kernelsize) - 0.5) * 2 * sqrt(6 / (fan_in + fan_out));
                    if p==1 
                        x=dbnmap(1:1,(25*j+1):(25*j+25));
                        net.layers{l}.k{i}{j}=reshape(x,5,5);  %(dbnmap(net.layers{l}.kernelsize^2*j-net.layers{l}.kernelsize^2+1:net.layers{l}.kernelsize^2*j)
                    else
                        x=dbnmap(150+(j-1)*6*49+49*(i-1)+1+25:150+(j-1)*6*49+49*i+25);
                        net.layers{l}.k{i}{j}=reshape(x,7,7); %(dbnmap(150+fan_in*(j-1)+1+net.layers{l}.kernelsize^2*i:150+fan_in*(j-1)+net.layers{l}.kernelsize^2*i+net.layers{l}.kernelsize)
        
                        end%卷积核的初始化生成一个7*7的卷积核，值为-1~1之间的随机数，再乘sqrt(6/(18*25))
                end
                net.layers{l}.b{j} = 0;
            end%卷积核初始化，1层卷积为1*6个卷积，2层卷积一共6*12=72个卷积核。
            inputmaps = net.layers{l}.outputmaps;
            p=2;
        end
   
    % 'onum' is the number of labels, that's why it is calculated using size(y, 1). If you have 20 labels so the output of the network will be 20 neurons.
    % 'fvnum' is the number of output neurons at the last layer, the layer just before the output layer.
    % 'ffb' is the biases of the output neurons.
    % 'ffW' is the weights between the last layer and the output neurons. Note that the last layer is fully connected to the output layer, that's why the size of the weights is (onum * fvnum)

    end
%%
   fvnum = prod(mapsize) * inputmaps;
    onum = size(y, 1);

    net.ffb = zeros(onum, 1);
    net.ffW = (rand(onum, fvnum) - 0.5) * 2 * sqrt(6 / (onum + fvnum));
end
