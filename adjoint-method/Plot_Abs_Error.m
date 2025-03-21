function Plot_Abs_Error(meshInfo,contourShow,iter_Evec,E_ext)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    nodsx = meshInfo.nelx+1;
    nodsy = meshInfo.nely+1;
    full_E = reshape(iter_Evec,nodsx,nodsy)';
    abs_error = 100*abs(E_ext-full_E)./E_ext;
    surf(contourShow.PlotX, contourShow.PlotY, abs_error); view(0,90);
    shading flat;
    shading interp;
    colormap('jet');
    colorbar;
    axis square;
    set(gca,'YDir','normal');
    axis off;
    pause(1);
end