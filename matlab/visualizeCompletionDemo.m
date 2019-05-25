%% draw trajectories
output_dir = '../output/2019-05-25-155717/';%../output/2019-05-25-155717/
dense_traj_name = [output_dir,  'dense_traj.txt'];
dense_traj= readTrajs(dense_traj_name) ;
  
drawTrajs(dense_traj, Inf,false,true,'uniform',200);
daspect([1,1,1])
set(gca,'XTick',[]);
set(gca,'YTick',[]);
hold off;
axis off;
 set(gcf,'position',[10,10,1000,1000])

addpath('plot2svg')
svgname = [output_dir, 'dense_traj.svg'];
plot2svg(svgname)
  