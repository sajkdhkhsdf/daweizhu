function P=algorithm1(M,N,P)
    ITER = 500;
    obj = zeros(ITER,1);
    [r,n]= size(P);
        %d = ones(r,1); %不用初始化
    d = 0.5./sqrt(sum((P.*P),2)+eps);%计算d，注意点乘，按行求和，不用转置了
for iter=1:ITER   
    D = spdiags(d,0,r,r);  %
    P = inv(D)*M'*inv(M*inv(D)*M')*N;
    d = 0.5./sqrt(sum((P.*P),2)+eps); 
    obj(iter) = sum(d);
    if (iter>1 && abs(obj(iter)-obj(iter-1))<0.00001)
        break;
    end
end;

