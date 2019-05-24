function [trajs,minX,minY,maxX,maxY] = readTraj(fname,N)
 if nargin < 2
     N=Inf
 end
 fname
f_id = fopen(fname,'r');
tline = fgetl(f_id);
 minX = Inf;minY =Inf;
 maxX = -Inf;maxY=-Inf;
%cluster_obj = {};
firstline = textscan(tline,'%d');
n=min( firstline{1},N);
fprintf('reading %d trajectories\n',n);
%n=5000;
trajs = cell(n,1);
for i=1:n
    tline = fgetl(f_id);
    if (length(tline)==1)
        continue;
    end 
    line = textscan(tline,'%s');
    tokens = arrayfun(@(x) str2num(char(x)),line{1});
     
    if tokens(1)== (length(tokens)-1)/2
        trajs{i} = reshape(tokens(2:end),2, tokens(1))';
    else
 %   trajs{i} = reshape(tokens(3:length(tokens)), 2,tokens(2) )';
     trajs{i} = reshape(tokens(3:length(tokens)), 2, (length(tokens)-2)/2  )';
    end
end

%drawTrajs(trajs,Inf)
end
