
function Q = f_GramSchmidt(W)

Q=zeros(size(W));

Q(:,1)=W(:,1)/norm(W(:,1));

for i=2:size(W,2)
   
    s_proj=0;
    for j=1:i-1
    s_proj=s_proj+((Q(:,j)'*W(:,i))/(Q(:,j)'*Q(:,j)))*Q(:,j);
    end
    
    e=W(:,i)-s_proj;
    Q(:,i)=e/norm(e);
end