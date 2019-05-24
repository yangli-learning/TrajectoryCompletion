###########################################################
# A script for handling output from trajectory completion #
###########################################################
source makeGraph.sh


mv -v dense_traj* ${dir}/.;
mv -v skeleton*.png ${dir}/.;
mv -v $jgraphDot ${dir}/.; 
mv -v completionResults.csv  ${dir}/.;
mv -v graph_clean*.dot ${dir}/.;
mv graph*.dot ${dir}/.;
cp -v ../GPStudio/parameters.ini ${dir}/.;
