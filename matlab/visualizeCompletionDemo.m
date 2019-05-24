%% draw trajectories
output_dir = '../output/2019-05-24-155615/';
%../output/2019-05-24-154241/';
%'../output/2019-05-24-153730/'; %good
%'../output/2019-05-24-153158/';
%'../output/2019-05-24-151559/';
%'../output/2019-05-24-151054/'
%'../output/2019-05-24-145102/'  ; %2019-05-24-140524
dense_traj_name = [output_dir, 'dense_traj.txt'];
 
dense_traj= readTraj(dense_traj_name) ;
 
%drawTrajsComp(dense_traj,[],Inf);
drawTrajs1(dense_traj, Inf,false,true,'uniform',150);
daspect([1,1,1])
set(gca,'XTick',[]);
set(gca,'YTick',[]);
hold off;
axis off;
 set(gcf,'position',[10,10,1000,1000])

addpath('plot2svg')
svgname = [output_dir, 'dense_traj.svg'];
plot2svg(svgname)

%% 
plot2svg('svg/region4A30.svg')


%% save to mat files 
na=[];
 
for i = 1:size(trajEdgeVote,1)
     na=char(trajEdgeVote(i,1));
     trajEdgeVote{i,3} = readTraj(na);
     fprintf('finish processing %s \n',na);
     
end;  
%%
 save('trajBest.mat','trajEdgeVote');

%% compare normalized density


figure(3);
ids= [1]%,2,3];
trajEdgeVoteEps=cell(length(ids),2);
for i =1:length(ids)
    id = ids(i);
   [trajEdgeVoteEps{i,1}, trajEdgeVoteEps{i,2}]=epsilonStatsRatio(trajEdgeVote{id,3});
end 
 
labels=cell(1,length(ids));
colors= linspecer(length(ids));  
for i=1:length(ids)
    id = ids(i);
    plot(trajEdgeVoteEps{i,2},trajEdgeVoteEps{i,1},'.-','MarkerSize',15,'LineWidth',1.5,'Color',colors(i,:));
    labels{i} =['Traj ', num2str(trajEdgeVote{id,2})]; 
        
    hold on;
end
legend(labels)

%% draw the best traj
bestId = 1;
N=20000;
drawTrajs( trajEdgeVote{bestId,3},N,false, true, 'a',100);

set(gcf,'color','w')
axis off
xlim([ 443000      450000])
ylim([4420500     4425000])
%% strip length
figure(3);
ids= [1,2,3];
trajEdgeVoteEps=cell(length(ids),2);
for i =1:length(ids)
    id = ids(i);
   [trajEdgeVoteEps{i,1}, trajEdgeVoteEps{i,2}]=epsilonStats(trajEdgeVote{id,3});
end 
 
labels=cell(1,length(ids));
colors= linspecer(length(ids));  
for i=1:length(ids)
    id = ids(i);
    plot(trajEdgeVoteEps{i,2},trajEdgeVoteEps{i,1},'.-','MarkerSize',15,'LineWidth',1.5,'Color',colors(i,:));
    labels{i} =['Rule ', num2str(trajEdgeVote{id,2})];
    %if mod(i,2)==1
    %    labels{i} = [num2str(trajEdgeVote{i,2}),'a'];
    %else
    %    labels{i} = [num2str(trajEdgeVote{i,2}),'b'];
    %end
        
    hold on;
end
legend(labels)
%%
%title('')
%title('Sampling distance of trajectories with N samples removed from endpoints');
xlabel('distance between consecutive trajectory points (\epsilon)')
ylabel('# of points with sampling points below \epsilon')
set(gca,'XTick',trajEdgeVoteEps{1,2})

set(gcf,'color','w')
%%
export_fig('compareEdgeVotes.pdf','-m2.5')

%%
%{ 
# discarded trajectories 21458 (0.351557 of 61037 all trajectories,0.810256 of 26483 non-empty trajectories
# discarded trajectories 16816 (0.275505 of 61037 all trajectories,0.796514 of 21112 non-empty trajectories
# discarded trajectories 16921 (0.277225 of 61037 all trajectories,0.794301 of 21303 non-empty trajectories
# of 100 dense trajecotires 
# discarded trajectories 23234 (0.380654 of 61037 all trajectories,0.877318 of 26483 non-empty trajectories
# discarded trajectories 18145 (0.297279 of 61037 all trajectories,0.859464 of 21112 non-empty trajectories
# discarded trajectories 18372 (0.300998 of 61037 all trajectories,0.862414 of 21303 non-empty trajectories>> 

%}