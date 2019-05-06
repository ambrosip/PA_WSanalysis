clc;
clear;
close all;

% =========================================================================
% ===== Descriptions ======================================================
% =========================================================================

Description = {...
	'This is a tutorial on how to use some of the more advanced features of distFig. Hit enter in the workspace to continue to the next step and generate five figures.'
	
	'Suppose you want to distribute the figures to the western part of the screen in two rows and two columns. Do so by executing:'
	'You can now only see figure #5 in the top left corner, which has been placed on top of figure #1. This is because the function restarts the pattern, when it runs out of space. The extra figure (#5) can be ignored by executing (We use the eastern part now):'
	'Note that you didn''t use the ''Only'' command, but instead you ignored the restart of the pattern. However, this can also be reproduced by executing either of the following, which do the same thing:'
	
	'Suppose you want to distribute a figure so it spans the middle two rows in a four row pattern. This cannot be done by the standard features, but can be done with the ''Adjust'' command. This feature moves the whole pattern by a specified ammount. Move all figure to the north-western part by executing:'
	'Now, move the whole pattern by 100 pixels to the right and 200 pixels down by executing:'
	'Instead of moving them a specific number of pixels, you can use the same command with scalars from 0 to 1, which are realtive sizes of the selected monitor size. Move the whole pattern half a screen down by executing:'
	'The function detect automatically whether the input are integers or scalars. Now - returning to the case where one figure was supposed to span the middle two rows of a four row pattern. This can thus be done by executing:'
	'And let us position figure #2 and #3 on top and below figure #1 by executing the following two commands:'
	'Note how the ''Offset'' feature was used to position figure #3 by skipping the two rows, where figure #1 was placed.'
	
	'distFig also features some minor function which are useful when positioning many figures (for instance 50). In such a case, it can be hard to see the actual figures due to the menu bar and the surrounding whitespace. The menu bars can be remove by executing:'
	'The whitespace can also be removed by making a tight fit around the figure. Do so by executing:'
	
	'Now there is a lot more figure and a lot less crap. :)'
	
	'The last feature of distFig is if you use Simulink. The figures from a Simulink model can be distributed by themself, along with the normal figures or completely ignored by executing one of the following:'
	'(The above functions were not executed now) - try them out afterwards to see how our five figures are not distributed when executing the first one, where only Simulink figures are distributed.'
	
	'Try and play around with the function by executing some of the used commands below:'
	'To see a description of the various features execute:'
	'Please contact me, if you find any bugs/errors on AndersSSimonsen@GMail.com or write a comment on File Exchange. :)'
	};

% =========================================================================
% ===== Commands ==========================================================
% =========================================================================

Command = {...
	{''}
	{'distFig(''Pos'',''W'',''Rows'',2,''Cols'',2'')'}
	{'distFig(''Pos'',''E'',''Rows'',2,''Cols'',2'',''Extra'',''Ignore'')'}
	{'distFig(''Pos'',''E'',''Rows'',2,''Cols'',2'',''Not'',5)','!distFig(''Pos'',''E'',''Rows'',2,''Cols'',2'',''Only'',[1,2,3,4])'}
	
	{'distFig(''Pos'',''NW'')'}
	{'distFig(''Pos'',''NW'',''Adjust'',[100,-200])'}
	{'distFig(''Pos'',''NW'',''Adjust'',[0,-0.5])'}
	{'distFig(''Pos'',''E'',''Adjust'',[0,-0.25],''Only'',1,''Rows'',2)'}
	{'distFig(''Pos'',''E'',''Only'',2,''Rows'',4)','distFig(''Pos'',''E'',''Only'',3,''Rows'',4,''Offset'',3)'}
	{''}
	
	{'distFig(''Pos'',''W'',''Menu'',''none'')'}
	{'distFig(''Pos'',''W'',''Tight'',true)'}
	{''}
	
	{'!distFig(''Pos'',''W'',''Simulink'',''only'')','!distFig(''Pos'',''W'',''Simulink'',''include'')','!distFig(''Pos'',''W'',''Simulink'',''ignore'')'}
	{''}
	
	{'>>>ALL<<<'}
	{'!help distFig'}
	{''}
	};

Trigger_All = numel(Command) - 2;
n = 0;
for i = 1:(Trigger_All - 1)
	for j = 1:numel(Command{i})
		n = n + 1;
		Command{Trigger_All}{n} = Command{i}{j};
	end
end

% =========================================================================
% ===== Steps =============================================================
% =========================================================================

% ===== TextWidth =========================================================
TextWidth = 60;

for i = 1:numel(Description)
	% ===== Description ===================================================
	fprintf('===== Step %02.0f/%02.0f %s\n',i,numel(Description),repmat('=',1,TextWidth - 17));
	Str = '';
	while (1)
		if (numel(Description{i}) >= TextWidth)
			NewLine = regexp(Description{i}(1:TextWidth),'[ .]');
			Str = strcat(Str,Description{i}(1:NewLine(end)),'\n');
			Description{i}(1:NewLine(end) - 0) = [];
		else
			Str = strcat(Str,Description{i},'\n');
			break;
		end
	end
	fprintf(Str);
	
	% ===== Command =======================================================
	if ~((isempty(Command{i}{1})) && (numel(Command{i}) == 1))
		for j = 1:numel(Command{i})
			if (numel(Command{i}{j}) ~= 0)
				if (~strcmp(Command{i}{j}(1),'!'))
					disp(sprintf('<a href="matlab:%s;">%s</a>',Command{i}{j},Command{i}{j})); %#ok<DSPS>
				else
					disp(sprintf('<a href="matlab:%s;">%s</a>',Command{i}{j}(2:end),Command{i}{j}(2:end))); %#ok<DSPS>
				end
			end
		end
	end
	
	% ===== Continue ======================================================
	input('');
	
	% ===== Execute =======================================================
	if (~((isempty(Command{i}{1})) && (numel(Command{i}) == 1)) && (i ~= Trigger_All))
		for j = 1:numel(Command{i})
			if (numel(Command{i}{j}) ~= 0)
				if (~strcmp(Command{i}{j}(1),'!'))
					eval(sprintf(Command{i}{j}));
				end
			end
		end
	end
	
	% ===== Generate figures ==============================================
	if (i == 1)
		for j = 1:5
			figure('Color','w');
			hold on;
			plot(linspace(-1,1,100),sin(linspace(-1,1,100) * pi * j) * 0.25 + 0.75,'-','LineWidth',2);
			plot(linspace(-1,1,100),-sin(linspace(-1,1,100) * pi * j) * 0.25 - 0.75,'-','LineWidth',2);
			text(0,0,num2str(j),'HorizontalAlignment','Center','FontSize',100,'Color',[0.9,0.9,1] * 0.25,'BackGroundColor','w');
			hold off;
			box on;
			xlim([-1,1]);
			ylim([-1,1]);
			xlabel('XLabel');
			ylabel('YLabel');
		end
	end
end