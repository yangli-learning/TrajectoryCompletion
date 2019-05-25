Trajectory Completion Tutorial
=========================
A tutorial on how to use GPStudio to complete trajectories

[TOC]



Input format
-------------------

### GPS trajectories

The input trajectories need to be stored in a protobuf (PBF) file. See  https://github.com/yangli1-stanford/gpsproto for details. 
In this tutorial we will use a set of artificial GPS trajectories located at `data/traj_pbf/region5C-large.trace-ss20-n5.pbf`.   

### Configuration file

Algorithm parameters are set in a configuration file `GPStudio/parameters.ini`. While it lists many parameters, in most cases only a few need to be tuned (which will be discuss later). 

We are also including a sample configuration file `GPStudio/parameters_beijing_taxi.ini` for some of the Beijing taxi data used in the paper. It should be adjusted depending on map size, street density, number of trajectories and the magnitude of GPS noise.

Output format
---------------------

The output dense trajectories are written as `.txt` files. The first line is the number of trajectories, 

	[number_of_trajectories] 

The rest of the file contains one trajectory per line, in the following format:

	[N] [x1] [y1] [x2] [y2] .... [xN] [yN]

where `N` is the number of points in this trajectory. `x` and `y` are positions in UTM coordinates. Several intermediate files are generated, including skeleton graphs and junction networks in `.dot` (graphviz) format. 

GUI Mode
--------------

The GUI allows interactive control of different steps of the algorithm. Inside the `/build` directory, type 
```
bin/GPStudio
```
You can adjust parameters on the fly using the Parameter editor (Hover the mouse over the **>>** sign on the right side of the toolbar, then click **Open Parameter Editor**. Press Ctrl-S or click **Update Paramters** to save changes after editing. 

To start, double click on the empty canvas and select **Open GPS Trajectories**. Open the sample PBF file `../data/traj_pbf/region5C-large.trace-ss20-n5.pbf`.

![opentrajectories](/home/yang/git/GPStudio-cleanup/doc/opentrajectories.png)


### Select Bounding Box
First we need to select the region of interest by dragging on the map while holding Ctrl.  If you need to start over,  delete the bounding box by Shift+Click on the box. 

Use mouse scrolling to zoom, hold the middle mouse button and drag to pan the view. 

Next click **Vote Map Image** on the menu bar to generate a grayscale *map image* of the selected region. 

![mapimage](/home/yang/git/GPStudio-cleanup/doc/mapimage.png)

Trajectory points and the map image can be toggled on and off using **Toggle Traj Points** and **Toggle Map Image**.

You can save the bounding box and its voting image by clicking **File > Save Bounding Box**. This will generate a `.yaml` file that can be loaded later.  

Useful parameters (defined in `GPStudio/parameters.ini`): 

* `map_image_resolution`: resolution of the grayscale map image  


### Extract skeleton 
First, click "initialize L1 Medial Skeleton" to generate skeleton points based on the density of the map image. 

![initialize](/home/yang/git/GPStudio-cleanup/doc/initialize.png)

Useful parameter:

- `skeleton_thresh`: increase this number to generate more skeleton points in low density regions.  

Clicking **>** once will perform one step of the L1 medial skeleton extraction iteration.   e.g. After 3-4 steps, points outside the junction will contract to a single line. 

![skeleton_extraction](/home/yang/git/GPStudio-cleanup/doc/skeleton_extraction.png)

Click **Extract Skeleton** to trace skeleton segments. 

![skeleton_branches](/home/yang/git/GPStudio-cleanup/doc/skeleton_branches.png)

For complex maps, you need to repeat **Extract Skeleton** multiple times to extract all skeleton branches.  Clicking on a skeleton point while holding Ctrl-Shift it will display the current support size of the selected point, and visualize the local variance in the lowerleft window.  You may also use **+** to increase the local support size (h) of each skeleton point if the skeleton points converge too slowly. 

Click **Compute Graph** to compute the final skeleton graph.

![skeleton_extraction](/home/yang/git/GPStudio-cleanup/doc/skeleton_graph.png)

Useful parameters:

* `skeleton_h0`: initial support size of L1-medial skeleton extraction. Increase this number to make skeleton extraction converge faster. 
* `skeleton_sigma_thresh`: decrease this value to include less salient skeletons branches

### Cluster junction trajectories 

Once the desirable skeleton graph is computed, use **Cluster junction trajectories** to compute flow clusters (trajectory segments grouped by traffic directions) . The detected clusters are shown inside the top left panel. In this example, a total of 12 clusters are discovered, with checkboxes that control the visibility of individual clusters. You can also filter the clusters by size using the slider. 

![flow_clusters](/home/yang/git/GPStudio-cleanup/doc/flowclusters.png)

Useful parameters:

`skeleton_proj_ratio_thresh`: a threshold that determines whether a projection is "stable". Increase this number to get more accurate trajectory flow clusters, but it will also reduce the number of completed trajectories.   

### Complete trajectories

After clusters are computed, click **Algorithm > Complete Trajectories** to generate dense trajectories. This process may take some time. Once completed, you can visualize the completion results of a specific trajectory by entering its id and click **Show Dense Trajectory**. For example, the image below shows the completed trajectory #400. Yellow arrows are the original GPS samples of the sparse trajectory. Red arrows are the added samples.

![completion](/home/yang/git/GPStudio-cleanup/doc/completion.png)

Command line mode
-----------------------------
The syntax is:
```
GPStudio  --interactive=false --trace_name [trajectory_file] --bbox [bounding_box_file]
```
The bounding box file can be generated using the GUI. A sample bounding box file is provided in `data/bbox/`. For example, 
```
bin/GPStudio  --interactive=false --trace_name ../data/traj_pbf/region5C-large.trace-ss20-n5.pbf --bbox ../data/bbox/region5C-large.trace-ss20-n5.bbox
```
The completed trajectory is named `dense_traj.txt`. You may also  use `build/guiless.sh` to run the program on command line. Its input is a list of trajectory and bounding box file names.  It will also move output files in a timestamped directory by the following naming convention: `output/[yyyy-mm-dd-hhmmss]/` . 

Visualize trajectories in MATLAB
-----------------------------------------------
A function `matlab/readTrac.m` is provided to read the output trajectories into MATLAB. 
To generate figures like Figure 1 in the paper, i.e. smoothed epsilon-dense trajectories, see `matlab/visualizeCompletionDemo.m`. 

