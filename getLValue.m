function bi = getLValue(i,tis,p,f)
[n,~] = size(tis);
bi = 0;

for m = 1:n
    tri = tis(m,:);
    bi = bi + abs(det(bsxfun(@minus,p(tri(2:end),:),p(tri(1),:))))/2;
end

bi = (1/6)*f(p(i,:))*bi;

end