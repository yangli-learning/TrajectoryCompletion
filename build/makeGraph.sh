###################################
# A script for visualizing graphs #
###################################

jgraphName=junction_graph_0
sgraphName=graph_clean
sgraphDot0="${sgraphName}0.dot"
sgraphDot1="${sgraphName}1.dot"
sgraphDot2="${sgraphName}2.dot"
sgraphDot3="${sgraphName}3.dot"
sgraphDot4="${sgraphName}4.dot"

jgraphDot="${jgraphName}.dot"

DATE=`date +%Y-%m-%d-%H%M%S`
dir="../output/$DATE"
if [ ! -d $dir ]; then
	mkdir -p $dir
fi
sgraphImg0="${dir}/${sgraphName}0.png" 
sgraphImg1="${dir}/${sgraphName}1.png" 
sgraphImg2="${dir}/${sgraphName}2.png" 
sgraphImg3="${dir}/${sgraphName}3.png" 
sgraphImg4="${dir}/${sgraphName}4.png" 

jgraphImg="${dir}/${jgraphName}.png"
echo "Writing $sgraphImg0"

dot -Kfdp -n -Tpng -Goverlap=true -o $sgraphImg0 $sgraphDot0 
echo "Writing $sgraphImg1" 

dot -Kfdp -n -Tpng -Goverlap=true -o $sgraphImg1 $sgraphDot1 

echo "Writing $sgraphImg2"

dot -Kfdp -n -Tpng -Goverlap=true -o $sgraphImg2 $sgraphDot2 

echo "Writing $sgraphImg3"

dot -Kfdp -n -Tpng -Goverlap=true -o $sgraphImg3 $sgraphDot3

echo "Writing $sgraphImg4"

dot -Kfdp -n -Tpng -Goverlap=true -o $sgraphImg4 $sgraphDot4

dot -Kfdp -n -Tpng  -Goverlap=true -o $jgraphImg $jgraphDot
echo "Wrote $jgraphImg"

#eog $jgraphImg
 
########################################################
#dotfiles=$(find . -name "*graph[0-9]*.dot")
#for f in $dotfiles; do
#	filename=$(basename $f)
#	pngName="${filename%.*}.png"
#	dname="./${filename%.*}O.dot"
    #gvpr  -o $dname "N[$.degree==0]{delete(0,$);}" $f

#	echo "read ${dname} write ${pngName}. " 
#	echo
#	eog -n ${pngName}  
#done  

