%X:d*n,n列数是样本个数，d是特征数
%Y：n*d
%W: d*c
function [W,b]=algorithm2(X,Y,lambda)
[d,n]=size(X);
[dd,c]=size(Y);
one=ones(n,1); %表示全1列向量
W=rand(d,c);
b=rand(c,1);
ITER=500;
obj = zeros(ITER,1);
 for iter=1:ITER  
     W1=X'*W+one*b'-Y; %n*c的矩阵
     d1 = 0.5./sqrt(sum((W1.*W1),2)+eps); %
     D1 = spdiags(d1,0,n,n); %D^ 与w1的行数相等，为n*n 
     I=eye(n); %单位矩阵
     s=1/(one'*D1*one);  %1/1'D^1
     Q = (I-s*(one*one')*D1')*(Y-X'*W);
     M = [(I-s*(one*one')*D1')*X',lambda*I]; 
     N = Y-s*(one*one')*D1'*Y; 
     P = [W;Q];
     P = algorithm1(M,N,P);
     W = P(1:d,:);
     b = s*(Y'-W'*X)*D1*one;
     obj(iter) = sum(b);
     if (iter>1 && abs(obj(iter)-obj(iter-1))<0.00001)
         break;
     end
end


% (W,b) = algorithm2(X,Y,100)







