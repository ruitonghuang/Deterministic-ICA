function sig = rescale(sig)


for i = 1:2
    mag = max(abs(sig(i,:)));
    sig(i,:) = sig(i,:)./mag;
    if sig(i,1) < 0
        sig(i,:) = -1.* sig(i,:);
    end
end
