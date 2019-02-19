function nn = dbnunfoldtonn(dbn, outputsize)
%DBNUNFOLDTONN Unfolds a DBN to a NN
%   dbnunfoldtonn(dbn, outputsize ) returns the unfolded dbn with a final
%   layer of size outputsize added.
    if(exist('outputsize','var'))
        %DBNUNFOLDTONN Unfolds a DBN to a NN  
        %outputsize�����Ŀ�����label��ǩ��������MINST����10��DBNֻ����ѧϰfeature  
        %����˵��ʼ��Weight����һ��unsupervised learning������supervised���ÿ�NN  
        size = [dbn.sizes outputsize];
    else
        size = [dbn.sizes];
    end
    nn = nnsetup(size); 
    %��ÿһ��չ�����Weight��ȥ��ʼ��NN��Weight  
    %ע��dbn.rbm{i}.c��ȥ��ʼ����bias���ֵ  
    for i = 1 : numel(dbn.rbm)
        nn.W{i} = [dbn.rbm{i}.c dbn.rbm{i}.W];
    end
end

