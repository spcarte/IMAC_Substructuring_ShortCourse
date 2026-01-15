function varargout = MAC_compare(phi1,phi2,varargin)
% Function to compute and plot the MAC between two sets of mode shapes.
% phi is assumed to be a matrix with mode shape vectors as columns.
%   varargin{1} = 2 or 3 to make a 2D or 3D plot of the MAC
%   varargin{2} = the figure number to use for the MAC plot
%   varargin{3} = frequencies of phi1
%   varargin{4} = frequencies of phi2

a = size(phi1,2); % Number of modes in Phi1
b = size(phi2,2); % Number of modes in Phi2
mac = zeros(a,b); % Preallocating MAC value matrix

% Loop over each mode combination
for ii=1:b
    for jj=1:a
        mac(jj,ii) = abs((phi1(:,jj)'*phi2(:,ii))^2 / ((phi1(:,jj)'*phi1(:,jj))*(phi2(:,ii)'*phi2(:,ii))));
    end
end

% Extra inputs define plot type
if ~isempty(varargin)
    if length(varargin)>1
        if isempty(varargin{2})
            figure
        else
            figure(varargin{2});
        end
        if length(varargin)>2
            lab1 = cellstr([num2str(round(varargin{3}(1:a))) repmat(' Hz',[a 1])]);
            lab2 = cellstr([num2str(round(varargin{4}(1:b))) repmat(' Hz',[b 1])]);
        else
            lab1 = cellstr(string(1:a)');
            lab2 = cellstr(string(1:b)');
        end
    end

    if varargin{1}==2
        % Plot the MAC Matrix as a 2D Image
        MACplot=pcolor(blkdiag(mac,0).');
        colormap(flipud(colormap('bone')));
        c=colorbar('eastoutside'); c.Label.String='MAC Value';
        xlabel('Phi 1'); ylabel('Phi 2'); title('MAC Plot');
        axis tight; axis equal;
        MACplot.Parent.XTick = (1:a)+.5;
        MACplot.Parent.XTickLabel = lab1;
        MACplot.Parent.YTick = (1:b)+.5;
        MACplot.Parent.YTickLabel = lab2;

    elseif varargin{1}==3
        % Plot the MAC Matrix as a 3D Bar Plot
        MACplot=bar3(mac,1);
        for ii = 1:length(MACplot)
            MACplot(ii).CData = MACplot(ii).ZData;
            MACplot(ii).FaceColor = 'interp';
        end
        colormap(flipud(colormap('bone')));
        c=colorbar('eastoutside'); c.Label.String='MAC Value';
        axis([.5 b+.5 .5 a+.5 -inf inf]); view([-90 90]);
        xlabel('Phi 2'); ylabel('Phi 1'); zlabel('MAC Value'); title('MAC Plot');
        MACplot(1).Parent.XTickLabel = lab2;
        MACplot(1).Parent.YTickLabel = lab1;
    else
        disp('Put 2 for 2D plot, 3 for 3D plot');
    end
end

% MAC matrix as output if asked for
if nargout==1
    varargout = {mac};
else
    varargout = {};
end

end
