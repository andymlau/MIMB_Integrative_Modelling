#~/imp-projects/imp-270314/imp-fast-mpi/setup_environment.sh

python='/anaconda2/bin/python'
initial_map=${1}
number_of_components=${2}

echo 'initial_map =' ${initial_map}
echo 'number_of_components = ' ${number_of_components}

# Format: python create_gmm.py map.mrc number_of_components output_textfile -s map_threshold -m out_map.mrc -a voxel_size -i number_of_GMM_iterations
# E.g. /anaconda2/bin/python create_gmm.py model_negstain_28032017.mrc 50 model_negstain_28032017.mrc.gmm.30.txt -s 0.05 -m model_negstain_28032017.mrc.gmm.30.mrc -a 4 -i 200
# ${python} ./create_gmm.py ${initial_map} ${2} ${initial_map%.mrc}.gmm.30.txt -s ${3} -m ${initial_map%.mrc}.gmm.30.mrc -a ${4} -i ${5}
${python} ./create_gmm.py ${initial_map} ${number_of_components} ${initial_map%.mrc}.gmm.n${number_of_components}.txt -m ${initial_map%.mrc}.gmm.n${number_of_components}.mrc
