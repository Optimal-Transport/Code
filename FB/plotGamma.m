function plotGamma(x,y,gamma,M,N)
% This function plots the transport plan gamma, where gamma_{i,j} is the
% mass transported from x_i to y_j

    figure
    image( [x(1),x(end)], [y(1),y(end)], gamma', 'CDataMapping', 'scaled' )
    set(gca, 'YDir', 'normal')
    colorbar
    axis equal
    axis tight
    xlabel('$x$', 'interpreter', 'latex', 'FontSize', 14)
    ylabel('$y$', 'interpreter', 'latex', 'FontSize', 14)
    title('Mass $\gamma_{i,j}$ transported from $x$ to $y$ with $M=' + string(M) + '$ $N=' + string(N) + '$', ...
        'interpreter', 'latex', 'FontSize', 14)

end