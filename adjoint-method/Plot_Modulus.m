function Plot_Modulus(meshInfo,contourShow,iter_Evec)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    nodsx = meshInfo.nelx+1;
    nodsy = meshInfo.nely+1;
    full_E = reshape(iter_Evec,nodsx,nodsy)';
    surf(contourShow.PlotX, contourShow.PlotY, full_E); view(0,90);
    shading flat;
    shading interp;
    colormap('jet');
    colorbar;
    axis square;
    set(gca,'YDir','normal');
    axis off;
    pause(1);
end