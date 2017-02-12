# to (re)-install

source $HOME/cmb/bin/python3_setup.sh

# for msprime
source /usr/usc/hdf5/1.8.12/setup.sh
export PKG_CONFIG_PATH=$HOME/cmb/msprime

python3 setup.py install --prefix /home/rcf-40/pralph/cmb/python3.5
