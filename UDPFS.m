function [W,A,WResult] = UDPFS(X,gamma,r,fea_num, rdim,c)
%UNTITLED5 ???????§Û?????????
%   ????????????
num = size(X,2);
dim = size(X,1);

X0 = X';
mX0 = mean(X0);
X1 = X0 - ones(num,1)*mX0;
scal = 1./sqrt(sum(X1.*X1)+eps);
scalMat = sparse(diag(scal));
X = X1*scalMat;
X=X';

M  = rand(dim,c); 
W  = rand(dim,rdim); 
A=rand(num,c);

for iter = 1:1
       %% complete clustering centroid
        F=W'*X;
        M=(A'*F')./(A'*ones(num,1));
        %% updating membership martix
        distx = L2_distance_1(F,M');

        for i=1:num
            dxi = distx(i,:);
            ad = -(dxi)/(2*r);
            A(i,:) = EProjSimplex_new(ad);
        end;
        %% complete with-in scatter matrix
         b=sum(A,2);
         dd=sum(A,1);
         dd=1./(dd+eps);
         B=diag(b);
         D=diag(dd);
         LL=B-A*D*A';
         Sw=X*LL*X';
         
         [W,WResult]=updateW(Sw,X,gamma,dim,rdim);%updateW(X, rdim, gamma, Sw);
end
end