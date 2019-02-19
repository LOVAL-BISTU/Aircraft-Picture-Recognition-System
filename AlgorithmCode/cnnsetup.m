function net = cnnsetup(net, x, y)
    assert(~isOctave() || compare_versions(OCTAVE_VERSION, '3.8.0', '>='), ['Octave 3.8.0 or greater is required for CNNs as there is a bug in convolution in previous versions. See http://savannah.gnu.org/bugs/?39314. Your version is ' myOctaveVersion]);
    inputmaps = 1;
    mapsize = size(squeeze(x(:, :, 1)));%��ȥsizeΪ1��ά��Ȼ���ȡ����������
    load dbnmap.mat%�޸ĺ����
    p=1;
    for l = 1 : numel(net.layers)   % layer��ÿһ�㡣numel�������ڼ�������������ָ��������Ԫ�ظ���������һ��ͼ����numel(A)������������������
        if strcmp(net.layers{l}.type, 's')%c = strcmp(str1,str2)���Ƚ��ַ��� str1 �� str2 ������ȫ����򷵻� 1 ������ȷ��� 0��
            mapsize = mapsize / net.layers{l}.scale;%���ְѲ���Ƶ�ʽ����������ǽ�������downsample�����������ĺô��Ǽ����������㣬Ҳ���Ǽ�������ʱ�䣬��ʵʱ����ʱ�����õķ������������scale=2������pooling֮��ͼ�Ĵ�С��pooling��֮��û���ص�������pooling���ͼ��Ϊ25*25.mapsize���������һ���size
            assert(all(floor(mapsize)==mapsize), ['Layer ' num2str(l) ' size must be integer. Actual: ' num2str(mapsize)]);
            for j = 1 : inputmaps
                net.layers{l}.b{j} = 0;
            end
        end
        if strcmp(net.layers{l}.type, 'c')%����Ǿ���㡣
            mapsize = mapsize - net.layers{l}.kernelsize + 1;%�õ�featuremap��size����СΪ�ϲ��С-����˴�С+1(��50-7+1=44����22-11+1����featruemap��ÿ���˲�������õ���ͬ�����ķ�ӳ��
            fan_out = net.layers{l}.outputmaps * net.layers{l}.kernelsize ^ 2;%�õ��ò��������ӵ�����=�������*�����size=6*��7*7����12��11*11��
            for j = 1 : net.layers{l}.outputmaps  %  output map �������*�����size
                fan_in = inputmaps * net.layers{l}.kernelsize ^ 2;     %�����һ�����map����Ӧ�����о���ˣ�������Ȩ����Ԫ������=1*7*7��6*11*11
                for i = 1 : inputmaps      %  input map%ԭCNNΪnet.layers{l}.k{i}{j} = (rand(net.layers{l}.kernelsize) - 0.5) * 2 * sqrt(6 / (fan_in + fan_out));
                    if p==1 
                        x=dbnmap(1:1,(25*j+1):(25*j+25));
                        net.layers{l}.k{i}{j}=reshape(x,5,5);  %(dbnmap(net.layers{l}.kernelsize^2*j-net.layers{l}.kernelsize^2+1:net.layers{l}.kernelsize^2*j)
                    else
                        x=dbnmap(150+(j-1)*6*49+49*(i-1)+1+25:150+(j-1)*6*49+49*i+25);
                        net.layers{l}.k{i}{j}=reshape(x,7,7); %(dbnmap(150+fan_in*(j-1)+1+net.layers{l}.kernelsize^2*i:150+fan_in*(j-1)+net.layers{l}.kernelsize^2*i+net.layers{l}.kernelsize)
        
                        end%����˵ĳ�ʼ������һ��7*7�ľ���ˣ�ֵΪ-1~1֮�����������ٳ�sqrt(6/(18*25))
                end
                net.layers{l}.b{j} = 0;
            end%����˳�ʼ����1����Ϊ1*6�������2����һ��6*12=72������ˡ�
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
