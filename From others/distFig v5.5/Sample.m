clc;
clear;
close all;

% =========================================================================
% ===== General figures ===================================================
% =========================================================================

Color = lines(7);
Color = Color([1,7,5,4],:);
Color(end+1,:) = [0.9,0.9,1] * 0.25; 
I = [1,2,2,2,3,3,3,4,4,4,4,5,5,5];
Size = [350,250,150,125,100];
for i = 1:numel(I)
	figure('Color','w');
	text(0,0,num2str(i),'HorizontalAlignment','Center','FontSize',Size(I(i)),'Color',Color(I(i),:));
	xlim([-1,1]);
	ylim([-1,1]);
	axis off;
	set(gcf,'Color',(Color(I(i),:) + [1,1,1]) / 2,'MenuBar','none');
	set(gca,'XTick',[],'YTick',[]);
end

% =========================================================================
% ===== Distribute figures ================================================
% =========================================================================

distFig('Pos','NW','Only',1,'Tight',true);
distFig('Pos','SW','Cols',3,'Only',(2:4),'Tight',true);
distFig('Pos','NE','Rows',2,'Only',5,'Tight',true);
distFig('Pos','NE','Rows',2,'Cols',2,'Only',[6,7],'Offset',2,'transpose',true,'Tight',true);
distFig('Pos','SE','Rows',2,'Cols',4,'Only',(8:11),'Tight',true);
distFig('Pos','SE','Rows',3,'Cols',2,'Not',(1:11),'Offset',3,'Tight',true);