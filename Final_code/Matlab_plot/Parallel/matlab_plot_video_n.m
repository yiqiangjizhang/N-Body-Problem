clc
clear all
close all

% Read csv file space separated
data = dlmread('output.csv');

%% Initialize video
myVideo = VideoWriter('myVideoFile'); %open video file
myVideo.FrameRate = 70;  %can adjust this, 5 - 10 works well for me
open(myVideo)

% Time
% t = data(:,1);
curve1 = animatedline('Marker','o','MaximumNumPoints',1,'Color','m','MarkerSize',3);
curve2 = animatedline('Marker','o','MaximumNumPoints', 1,'Color','r','LineWidth',3);
curve3 = animatedline('Marker','o','MaximumNumPoints', 1,'Color','r','LineWidth',3);
curve4 = animatedline('Marker','o','MaximumNumPoints', 1,'Color','b','LineWidth',3);
curve5 = animatedline('Marker','o','MaximumNumPoints', 1,'Color','r','LineWidth',3);
curve6 = animatedline('Marker','o','MaximumNumPoints', 1,'Color','g','LineWidth',3);
curve7 = animatedline('Marker','o','MaximumNumPoints', 1,'Color','k','LineWidth',3);
curve8 = animatedline('Marker','o','MaximumNumPoints', 1,'Color','c','LineWidth',3);
curve9 = animatedline('Marker','o','MaximumNumPoints', 1,'Color','b','LineWidth',3);
curve10 = animatedline('Marker','o','MaximumNumPoints', 1,'Color','k','LineWidth',3);
curve = [curve1, curve2, curve3, curve4, curve5, curve6, curve7, curve8, curve9, curve10];


set(gca, 'XLim', [-6e9 6e9], 'YLim', [-6e9 6e9]);
box on;
grid on;
grid minor;

for j = 1:500:size(data,1)
    counter = 1;
    
    for i=2:6:size(data,2)
        
        addpoints(curve(counter),data(j,i),data(j,i+1));
        drawnow;
        hold on;
        counter = counter + 1;
        frame = getframe(gcf); %get frame
        writeVideo(myVideo, frame);
    end
    str = {strcat('Days = ', num2str(round(data(j,1)/86400)))};
    legend(str{:});
end
close(myVideo);

