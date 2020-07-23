function [W,WResult] = updateW(Sw,X,gamma,dim,m)
% function  [W,obj]=InterationW(X, c, gamma,Sw)
%L: laplacian matrix
%X: data matrix(dim*num)
%gamma: coefficient of L21
%dim: dimension of X
%m: projection dimension of W

INTER_W = 2;
Q = eye(dim);
xlx= Sw;
xlx = (xlx+xlx')/2;
p=1; % L_2p
for i = 1:INTER_W
    tempXLXQ = (xlx+gamma*Q);
    [vec,val] = eig(tempXLXQ);
    [~,di] = sort(diag(val));
    W = vec(:,di(1:m));

    tempQ = 0.5*p * (sqrt(sum(W.^2,2)+eps)).^(p-2);
    Q = diag(tempQ);

    w1(i) = trace(W'*Sw*W); % log Tr(WXLXW)
    w2(i) = gamma*sum(sqrt(sum(W.^2,2)));% gama*||W||_21
    WResult(i) = w1(i)+w2(i);

    if i > 1 && abs(WResult(i-1)-WResult(i)) < 0.000001
        break;
    end;
    
end;
% WResult = WResult';
% plot(WResult)
end