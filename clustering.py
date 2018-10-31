import sys
#sys.path.append('/opt/imp-2.7.0/lib')
sys.path.append('/data/users/shared/software/anaconda/anaconda2/lib')
sys.path.append('/data/users/shared/software/anaconda/anaconda2/lib/python2.7/site-packages')

import matplotlib as mpl
mpl.use('Agg')

import IMP
import IMP.pmi
import IMP.pmi.macros
import os

# most common settings
num_clusters = 5
num_top_models = 100
merge_directories = ["./"]
#"/Users/andylau/Dropbox/COP9-CRL2/CRL2_XLMS/Nedd8_rotational_transition/IMP/IMP_N1_CRL2~N8_4P5O_pre-neddylation_positions/"
prefiltervalue = 2900.0
out_dir = "kmeans_%i_%i/" %(num_top_models,num_clusters)
if '--test' in sys.argv: prefiltervalue=8000.0

#################################
# should not have to change below
##################################

model=IMP.Model()

# initialize the macro
mc=IMP.pmi.macros.AnalysisReplicaExchange0(model,
                                           merge_directories=merge_directories)

# fields that have to be extracted for the stat file
feature_list=["ISDCrossLinkMS_Distance_intrarb",
              "ISDCrossLinkMS_Distance_interrb",
              "ISDCrossLinkMS_Data_Score",
              "GaussianEMRestraint_None",
              "SimplifiedModel_Linker_Score_None",
              "ISDCrossLinkMS_Psi",
              "ISDCrossLinkMS_Sigma"]

# Dictionary of densities to be calculated
              
density_names = {"TS_A":["TS_A"],"TS_B":["TS_B"],"TS_C":["TS_C"],"TS_D":["TS_D"]}

# list of component names needed to calculate the RMSD for the clustering
rmsd_names = {"TS_A":"TS_A",
              "TS_B":"TS_B",
              "TS_C":"TS_C",
              "TS_D":"TS_D"}

# components used for structural alignment
align_names = None # (None because EM provides reference frame)

mc.clustering(prefiltervalue=prefiltervalue,                   # prefilter the models by score
              number_of_best_scoring_models=num_top_models,    # number of models to be clustered
              alignment_components=None,                       # list of proteins you want to use for structural alignment
              rmsd_calculation_components=rmsd_names,          # list of proteins used to calculated the rmsd
              distance_matrix_file="distance.rawmatrix.pkl",   # save the distance matrix
              outputdir=out_dir,                               # location for clustering results
              feature_keys=feature_list,                       # extract these fields from the stat file
              load_distance_matrix_file=False,                 # skip the matrix calculation and read the precalculated matrix
              display_plot=True,                               # display the heat map plot of the distance matrix
              exit_after_display=False,                        # exit after having displayed the distance matrix plot
              get_every=1,                                     # skip structures for faster computation
              number_of_clusters=num_clusters)             # setup the list of densities to be calculated
