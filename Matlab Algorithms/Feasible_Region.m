x = linspace(0.9,1.05);
y_2 = x + 1 - (3/4)^0.5;
y_4 = x + (3/4)^0.5 - (1/2)^0.5;


plot(x,y_2)
hold on
plot(x,y_4)
yline((5/4)^0.5);
yline(1-(5/4)^0.5+(3/2)^0.5);
xline((1/2)^0.5-1+(3/2)^0.5);
hold off