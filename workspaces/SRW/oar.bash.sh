#!/bin/bash -l
#OAR -l nodes=1/core=28,walltime=4

CMD_LAT="python config_03_SE.py "$@

cat $OAR_NODEFILE
which mpirun

echo PATH=$PATH
env | grep PYTHON

# INFINIBAND: -q ib    (hosts like hib*)
# CMD=$(echo mpirun -machinefile $OAR_NODEFILE -n $(wc -l $OAR_NODEFILE | awk '{print $1}') \
#      --mca btl openib,self --mca btl_tcp_if_exclude ib0,lo,eth0,eth1,docker0 \
#      --mca orte_rsh_agent  "oarsh" $CMD_LAT)

CMD=$(echo mpiexec -machinefile $OAR_NODEFILE -n $(wc -l $OAR_NODEFILE | awk '{print $1}') \
      -env btl tcp,self -env btl_tcp_if_exclude ib0,lo,docker0 \
      -launcher ssh -launcher-exec /usr/bin/oarsh $CMD_LAT )

echo $CMD
$CMD
