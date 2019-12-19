load exampledata


opts.k = 20;
opts.branch = 3;
opts.tip = 0.2;

[I1,I2,I3] = findpoint(X,opts);

subplot(1,2,1)
plot(X(1,:),X(2,:),'m.'),hold on

subplot(1,2,2)
plot(X(1,I1),X(2,I1),'g.',...
            X(1,I2),X(2,I2),'k.',...
            X(1,I3),X(2,I3),'r.')