function errorbarplot_sb(X,Y,color1,color2,opt)
if nargin<3
    color1='k';
end
if nargin<4
    color2=[.5 .5 .5];
end
% if nargin<5
%     xlims=[-5 5];
% end
X=X(:)';
% Ymean=smooth(nanmean(Y,1),5)';
if nargin<5 || isempty(opt)
    Ymean=nanmean(Y,1);
    y1=Ymean+(nanstd(Y,[],1))./sqrt(sum(~isnan(Y),1));
    y2=Ymean-(nanstd(Y,[],1))./sqrt(sum(~isnan(Y),1));
    
elseif opt==1
    Ymean=nanmedian(Y,1);    
    y1=Ymean+(nanstd(Y,[],1));
    y2=Ymean-(nanstd(Y,[],1));
  
else
    y1=(nanmax(Y,[],1));
    y2=(nanmin(Y,[],1));
end
    X2 = [X fliplr(X)];
    Y = [y1 fliplr(y2)];
    hold on
    fill(X2,Y,color1,'LineStyle','none','LineWidth',1,'FaceAlpha',.2);
    hold on
    plot(X,Ymean,'-','Color',color2)
% for rep=1:length(X)
%     plot([X(rep) X(rep)],[y1(rep) y2(rep)],'k-')
% end
% plot(X,Ymean,'k-')