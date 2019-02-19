function [result,net]=cnnrecg(test_x,tag)
if tag==1
    load mynet2;
    net=cnnff(net,test_x);
    a=net.o;
    a=round(a);
    result=a(1); 
else
    load fenlei;
    net=cnnff(net,test_x);
    b=net.o;
    b=round(b);
    result=b(1);
end
end

