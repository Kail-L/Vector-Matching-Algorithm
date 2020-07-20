function [R]=gschmidt(Phi,T,W,P)
[m,n]=size(P);
R(1,:) = Phi'*T*sqrtm(W)./sqrt(Phi'*P*Phi);
e = eye(n);
for i=2:n
    ek = e(:,i);
    S = 0;
    for j=1:i-1
        S = S+(ek'*R(j,:)')*R(j,:)';
    end
    R(i,:) = ek - S;
    if R(i,:) == 0
        R(i,:)= e(:,1) - S;
    end
    R(i,:) = R(i,:)/norm(R(i,:));
end
end