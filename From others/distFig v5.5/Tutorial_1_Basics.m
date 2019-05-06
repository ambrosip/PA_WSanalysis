clc;
clear;
close all;

% =========================================================================
% ===== Descriptions ======================================================
% =========================================================================

Description = {...
	'This is a tutorial on how to use distFig. Hit enter in the workspace to continue to the next step and generate five figures.'
	
	'By executing distFig(); with no arguments, the function will distribute all open and visible figures on the central screen. Hit enter to execute:'
	'Various arguments can be passed to the function to distribute the figures to different parts of the screen(s). Hit enter to distribute the figures in the western part of the screen by executing:'
	'Suppose you want to distribute the figures in five rows instead. Hit enter to execute:'
	'You can also distribute them in five columns instead. Hit enter to execute:'
	'You can also choose which figures to distribute. Hit enter to only distribute figure #1, #3 and #4 by executing:'
	'These can also be distributed in three rows by adding (''Rows'',3) to the function. Hit enter to execute:'
	'You can also choose to not distribute specific figures. Hit enter to distribute all but figure #5 by executing:'
	'Note how the figures are distributed in the columns at first, after which they restart in the next column in the top. Use (''Transpose'',true) to transpose the order, at which they are distributed. Hit enter to execute:'
	
	'Suppose you only want to distriubute figure #1 to a specific place in the 2x2 pattern. Use input argument ''Offset'' to do so. We will use the eastern part of the screen for this. Hit enter to execute:'
	'Hit enter to offset it by 1 by executing:'
	'Hit enter to offset it by 2 by executing:'
	'Hit enter to offset it by 3 by executing:'
	'In order to do so, it is recommended to supply the number of column or rows, as the pattern will then be fixed. Otherwise the function will optimize the aspect ratio. This can be demonstrated by executing:'
	'Now figure #1 fills the whole area (''W''), but if all figures are distributed to this part of the screen, the pattern will change.  Hit enter to execute:'
	'To demonstate the ''Offset'' feature again we add (''Offset'',1). Hit enter to execute:'
	
	'Suppose you want to distribute figure #1 to the left and the remainding to the right. Do so by executing distFig two times using the ''Only'', ''Not'' and ''Offset'' arguments. Hit enter to execute:'
	'And the other way around by executing:'
	
	'Get it? It might seem a little complicated, but it''ll catch on to you. If you have an external monitor you can distribute it by using ''Screen''. Hit enter to execute:'
	'You might have a warning now, if you don''t have an external monitor. You can also pass the direction of the monitor, if you have multiple external monitors.'
	
	'Try and play around with the function by executing some of the used commands below:'
	'The function has other features as well, which are described in the function. To see a description of the various features execute:'
	'Please contact me, if you find any bugs/errors on AndersSSimonsen@GMail.com or write a comment on File Exchange. :)'
	};

% =========================================================================
% ===== Commands ==========================================================
% =========================================================================

Command = {...
	{''}
	
	{'distFig();'}
	{'distFig(''Pos'',''W'');'}
	{'distFig(''Pos'',''W'',''Rows'',5);'}
	{'distFig(''Pos'',''W'',''Columns'',5);'}
	{'distFig(''Pos'',''W'',''Only'',[1,3,4]);'}
	{'distFig(''Pos'',''W'',''Only'',[1,3,4],''Rows'',3);'}
	{'distFig(''Pos'',''W'',''Not'',5);'}
	{'distFig(''Pos'',''W'',''Not'',5,''Transpose'',true);'}
	
	{'distFig(''Pos'',''E'',''Only'',1,''Offset'',0,''Rows'',2,''Columns'',2);'}
	{'distFig(''Pos'',''E'',''Only'',1,''Offset'',1,''Rows'',2,''Columns'',2);'}
	{'distFig(''Pos'',''E'',''Only'',1,''Offset'',2,''Rows'',2,''Columns'',2);'}
	{'distFig(''Pos'',''E'',''Only'',1,''Offset'',3,''Rows'',2,''Columns'',2);'}
	{'distFig(''Pos'',''W'',''Only'',1);'}
	{'distFig(''Pos'',''W'');'}
	{'distFig(''Pos'',''W'',''Offset'',1);'}
	
	{'distFig(''Pos'',''W'',''Only'',1,''Columns'',2);','distFig(''Pos'',''W'',''Not'',1,''Rows'',4,''Columns'',2,''Offset'',4);'}
	{'distFig(''Pos'',''W'',''Only'',1,''Rows'',2);','distFig(''Pos'',''W'',''Not'',1,''Rows'',2,''Columns'',4,''Offset'',4,''Transpose'',true);'}
	
	{'distFig(''Screen'',''External'');'}
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
			text(0,0,num2str(j),'HorizontalAlignment','Center','FontSize',100,'Color',[0.9,0.9,1] * 0.25);
			xlim([-1,1]);
			ylim([-1,1]);
			axis off;
			set(gca,'XTick',[],'YTick',[]);
		end
	end
end