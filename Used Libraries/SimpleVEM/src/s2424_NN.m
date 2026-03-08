function [NN] = s2424_NN(pts,evaluatePts)

p1 = pts(1,:);
p2 = pts(2,:);
p3 = pts(3,:);
p4 = pts(4,:);

n = [ 1 p1(1) p1(2) p1(1)*p1(2) ; ...
    1 p2(1) p2(2) p2(1)*p2(2) ; ...
    1 p3(1) p3(2) p3(1)*p3(2) ; ...
    1 p4(1) p4(2) p4(1)*p4(2) ];

NN = zeros(1,length(n));

x = evaluatePts(1);
y = evaluatePts(2);

if rank(n) < length(n)
    % FORCE THE SLUTION WITH pinv():
    % N = pinv(n);
    % NN(i,:)=[1 x y x*y]*N;
    disp('[ERROR] Non compatible element!');
else
    N = inv(n);
    NN=[1 x y x*y]/n;
end

end