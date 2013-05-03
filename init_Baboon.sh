##########################################
# init_Baboon.sh file to load environment
# @author : rete , IPNL / CNRS
# @date : 03/05/13
#!/bin/sh
##########################################


export BABOON_HOME="/gridgroup/ilc/rete/Baboon"

export BABOON_BIN_DIR=$BABOON_HOME"/bin"
export BABOON_LIB_DIR=$BABOON_HOME"/lib"

export PATH=$PATH:$BABOON_BIN_DIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BABOON_LIB_DIR

