function [I1,I2,I3] = findpoint(X,opts)
%parameter introduction
% X   m*n gene expression matrix ,m genes and n cell 
% opts.k   the number of the neighbour 

k = opts.k;
branch = opts.branch;
tip = opts.tip;

[m,n] = size(X);
label = zeros(n,1);

Dis = DISMAT(X');
[D,I] = sort(Dis);

for i = 1: size(X,2)
    Ni = X(:,I(1:k,i));
    centeri = mean(Ni,2);
    CNi = Ni - repmat(centeri,1,k);
    [Ui,Si,Vi] = svd(CNi);
    sigma(:,i) = diag(Si);
    d(i)=sigma(1,i)/sigma(2,i);
end

I3 = find(d<branch);
I12 = find(d>=branch);

for i = 1:size(I12,2)
    I12(i);
    Ni = X(:,I(1:k,I12(i)));
    for j = 2:k
        vector(:,j) = Ni(:,j)-X(:,I12(:,i));
    end
    for m = 2:k
        for n = 2:k
            cos{i}(m,n)= vector(:,m)'*vector(:,n)/norm(vector(:,m))/norm(vector(:,n));
        end
    end
end

for i = 1:size(I12,2)
    cosi = cos{i};
    negative(i) = (size(find(cosi<0),1))/2;
    positive(i) = (k-2)*(k-1)-negative(i);
    rate(i) = negative(i)/((k-1)*(k-2)/2);
end

I1 = find(rate<tip);
I1 = I12(I1);
label(I3) = 3;
label(I1) = 1;
I2 = find(label==0);
label(I2) = 2;

end

