load('figure_difass.mat')
x= 15:5:100;
y=4:21;

figure(1)
plot(x,drawPointMtx1(2,y),'-+',...
    x,drawPointMtx1(3,y),'-o',...
    x,drawPointMtx1(4,y),'-x',...
    x,drawPointMtx1(5,y),'-s',...
   'LineWidth',2,'MarkerSize',8);
legend('LTE','Proposed-200pilot','Proposed-200pilot-PPE','Proposed-200pilot-PPCE');
xlabel('Number of active users')
ylabel('Pilot detection correct ratio')
grid on


