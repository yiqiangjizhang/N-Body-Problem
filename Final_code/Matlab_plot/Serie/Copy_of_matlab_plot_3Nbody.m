clc
clear all
close all

% Read csv file space separated
data= dlmread('output0.csv');
data2=dlmread('output1.csv');
data3=dlmread('output2.csv');

%% Initialize video
myVideo = VideoWriter('myVideoFile'); %open video file
myVideo.FrameRate = 5;  %can adjust this, 5 - 10 works well for me
open(myVideo)


% Time
% t = data(:,1);
curve1 = animatedline('Marker','o','MaximumNumPoints',1,'Color','m','MarkerSize',3);
curve2 = animatedline('Marker','o','MaximumNumPoints', 1,'Color','r','LineWidth',3);
curve3 = animatedline('Marker','o','MaximumNumPoints', 1,'Color','r','LineWidth',3);

curve = [curve1];
curve_2 = [curve2];
curve_3 = [curve3];

set(gca, 'XLim', [-15e7 15e7], 'YLim', [-15e7 15e7]);
% set(gca, 'XLim', [-6e10 6e10], 'YLim', [-6e10 6e10]);
box on;
grid on;
grid minor;


% For each time step
for j = 1:1:size(data,1)
    counter = 1;
    % For each planet
    for i=2:6:size(data,2)
        addpoints(curve(counter),data(j,i),data(j,i+1));
        drawnow;
        hold on;
        addpoints(curve_2(counter),data2(j,i),data2(j,i+1));
        drawnow;
        hold on;
        addpoints(curve_3(counter),data3(j,i),data3(j,i+1));
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