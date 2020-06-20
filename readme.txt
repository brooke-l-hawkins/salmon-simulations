
----------------------------------------------------------------------------------
Individual-based simulations suggest mixed impacts of warmer temperatures and a non-native predator on Chinook salmon

B.L. Hawkins, A.H. Fullerton, B.L. Sanderson, and E.A. Steel
----------------------------------------------------------------------------------
STEP I: Set up.

1) Download R and RStudio.

2) Download project files from GitHub at https://github.com/brooke-l-hawkins/salmon-simulations.

3) Unzip project files within salmon-simulations.zip.

4) Create new directory within salmon-simulations called 'data.out'.

5) Check that your files are stored with the following directory structure within your R project directory.

   code
      figure7.R
      figures4-6.R
      functions.R
      model.R
      summary.R
      sensitivity_analysis.R

   data.in
      fish.growth.lookup
         wt.growth.array.chinook.RData
         wt.growth.array.lmb.RData
      parameters
         parameters.csv
      shapefiles
         Basin_snq.dbf
         Basin_snq.sbx
         Basin_snq.shx
         Basin_snq.prj
         Basin_snq.shp
         Basin_snq.sbn
         Basin_snq.shp.xml
      sno.ssn
         WT.df2014.csv
         WT.df2015.csv
         binaryID.db
         chin_spawn_1120.dbf
         chin_spawn_1120.prj
         chin_spawn_1120.shp
         chin_spawn_1120.shx
         chin_spawn_1792.dbf
         chin_spawn_1792.prj
         chin_spawn_1792.shp
         chin_spawn_1792.shx
         chin_spawn_2241.dbf
         chin_spawn_2241.prj
         chin_spawn_2241.shp
         chin_spawn_2241.shx
         chin_spawn_2689.dbf
         chin_spawn_2689.prj
         chin_spawn_2689.shp
         chin_spawn_2689.shx
         chin_spawn_3361.dbf
         chin_spawn_3361.prj
         chin_spawn_3361.shp
         chin_spawn_3361.shx
         distance
            obs
               dist.net1.RData
            preds
               dist.net1.RData
               dist.net1.a.RData
               dist.net1.b.RData
         dnsegs.RData
         edges.cpg
         edges.dbf
         edges.prj
         edges.shp
         edges.shp.xml
         edges.shx
         jct.list.RData
         lmb_prevalent_250.dbf
         lmb_prevalent_250.prj
         lmb_prevalent_250.shp
         lmb_prevalent_250.shx
         lmb_prevalent_400.dbf
         lmb_prevalent_400.prj
         lmb_prevalent_400.shp
         lmb_prevalent_400.shx
         lmb_prevalent_500.dbf
         lmb_prevalent_500.prj
         lmb_prevalent_500.shp
         lmb_prevalent_500.shx
         lmb_prevalent_600.dbf
         lmb_prevalent_600.prj
         lmb_prevalent_600.shp
         lmb_prevalent_600.shx
         lmb_prevalent_750.dbf
         lmb_prevalent_750.prj
         lmb_prevalent_750.shp
         lmb_prevalent_750.shx
         netID1.dat
         preds.cpg
         preds.dbf
         preds.prj
         preds.sbn
         preds.sbx
         preds.shp
         preds.shp.xml
         preds.shx
         sites.cpg
         sites.dbf
         sites.prj
         sites.sbn
         sites.sbx
         sites.shp
         sites.shp.xml
         sites.shx
         upsegs.RData

   data.out

6) Open salmon-simulations R project in RStudio.

7) Install libraries in the setup section of model.R.
----------------------------------------------------------------------------------
STEP II: Run simulations for four scenarios.

To run each simulation, you will need to update 'climate.year' and 'predation.flag' variables in the model.R code, and then run the script.

1) Run Baseline scenario
      climate.year = 2014
      predation.flag = FALSE

2) Run Warm scenario
      climate.year = 2015
      predation.flag = FALSE

3) Run Predator scenario
      climate.year = 2014
      predation.flag = TRUE

4) Run Warm-Predator scenario
      climate.year = 2015
      predation.flag = TRUE

After running all four scenarios, you should have new directories and files in your data.out folder.

Each sub-directory corresponds to one scenario. 'A' indicates predators were absent, 'P' indicates predators were present; '2014' and '2015' indicate the climate year.

   data.out
      A.2014
         images
         quick.summary.2014.1.txt
         run.info.2014.1.txt
         salmon.array.2014.1.RData
         salmon.finalstep.2014.1.RData
      A.2015
         images
         quick.summary.2015.1.txt
         run.info.2015.1.txt
         salmon.array.2015.1.RData
         salmon.finalstep.2015.1.RData
      P.2014
         fish_other.array.2014.1.RData
         fish_other.finalstep.2014.1.RData
         images
         quick.summary.2014.1.txt
         run.info.2014.1.txt
         salmon.array.2014.1.RData
         salmon.finalstep.2014.1.RData
      P.2015
         fish_other.array.2015.1.RData
         fish_other.finalstep.2015.1.RData
         images
         quick.summary.2015.1.txt
         run.info.2015.1.txt
         salmon.array.2015.1.RData
         salmon.finalstep.2015.1.RData

The 'images' directories each contain 730 PNG files that picture the Snoqualmie basin and simulated fish. Black points represent where Chinook salmon will be spawned, red points represent Chinook salmon that have been spawned and are incubating, and blue points represent Chinook salmon that have emerged and are moving around the watershed. Gray boxes represent largemouth bass. Points disappear when fish die or outmigrate.

5) Create figures 4-7.

Run script 'figures4-6.R'. Run script 'figure7.R'.

After running these two scripts, you should have one new directory and eight new files in your data.out folder.

   data.out
      Figures
         Figure4.png
         Figure5Subyearling.png
         Figure5Yearling.png
         Figure6.png
         Figure7Spatial.png
         Figure7Temporal.png

Figure 5 is a combination of Figure5Subyearling.png and Figure5Yearling. Figure 7 is a combination of Figure7Spatial.png and Figure7Temporal.png.

----------------------------------------------------------------------------------

STEP III: Run global sensitivity analysis.

1) Install libraries in the setup section of sensitivity_analysis.R.

To run simulations for sensitivity analysis, run the ‘model.R’ script with the global variable 'SA' set to TRUE and the global variable ‘iter.list’ set to the number of iterations you wish to run (i.e., 1:300). The relevant function employed to vary parameters simultaneously is ‘fncGetParametersSA’ in the ‘functions.R’ file.

2) Run sensitivity analysis
      SA = TRUE
      iter.list = 1:300

To speed up the process, you may wish to break simulations up into smaller portions (e.g., iter.list = 1:25, 26:50, etc.) in order to simultaneously run multiple instances of R on one or more computers.

Once all simulations are completed, the results can be collected, analyzed, and visualized using random forest with the script ‘sensitivity_analysis.R’. There, you can:
 
     a) Collate parameters sets into one matrix
     b) Collate simulated data output metrics
     c) Use random forest regression to evaluate parameter sensitivity
     d) Produce plots to visualize sensitivity analysis results, and
     e) Replicate Figure 8 and all figures in Appendix 2

3) Run the sensitivty_analysis.R script line by line. This script is best used by walking through incrementally instead of sourcing, since you must choose which scenarios to interact with, and because some steps can take a long time to process and may be deemed unnecessary.


