

function [epsilon_stats,epsilons] = drawTrajs(trajs,N,randOn, smoothOn, ...
    colorscheme,minGap)
   figure();
 
 
  epsilons=0:30:150;
   epsilon_stats=zeros(1,length(epsilons));
  if nargin < 6
      minGap = 30000;
      if nargin < 5
          colorscheme = 'uniform';
          if nargin < 4
              smoothOn=false;
              if nargin < 3
                  randOn=false;
              end
          end
      end
  end
  
  n = min(N,length(trajs));
  if (randOn)
      colors = jet(n+1);
      idx =   randi(length(trajs),n,1)
  else  
      
      colors=jet(ceil(max(300,minGap)));
      idx =  1:1:n;
  end 
    nError =0;
      set(gcf,'color','w');
    for i=1:n 
        t = idx(i);
        raw = trajs{t};
        if (size(raw,1)<=5)
            continue;
        end
        xs = raw(:,2);ys=raw(:,1);
        if (~randOn)
            dx = diff(xs); dy = diff(ys);
            diffs = sqrt(dx.^2 +dy.^2);
            dd =max(diffs);
            epsilon_stats = epsilon_stats + cumsum(hist(diffs,epsilons));
            if( dd>minGap   )
                nError=nError + 1;
                continue;
            end
        end
        if (smoothOn)
           index = [1:1:length(xs)];
           xs = smooth(index,xs,12, 'lowess'); 
           ys = smooth(index,ys,12,'lowess'); 
         
        
        end;
           if (randOn)
               color = colors(1+mod(idx(i),size(colors,1)) ,:);
           else
               color=colors(max(1,round(dd)),:);
           end;
         if strcmp( colorscheme,'uniform')
            plot(ys,xs,'-','Color',[0    0.4470    0.7410], 'LineWidth',0.1);
         else
            plot(ys,xs,'-','Color',color,'LineWidth',0.1);
         end
         hold on;
    end
    fprintf('# discarded trajectories %d (%f of %d trajectories)', ...
    nError,nError/n,n); 
      set(gcf,'color','w');
end