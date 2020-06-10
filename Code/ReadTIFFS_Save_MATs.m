clear

dir0 = dir('*_Timelapse_???');

numDirs = size(dir0,1);
%%
for k=5 %1:numDirs   %15:16%7:12%numDirs
    disp(dir0(k).name)
    readNeutrophils(dir0(k).name);
end
disp('------------ FINISHED ------------')

beep

%%
%k=5;
%disp(dir0(k).name)
%readNeutrophils(dir0(k).name);