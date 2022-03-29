#!/bin/bash

echo 'Setting up simulation core' $1 'with simulation id' $2;
rm -v MCAS_v2_sim$1/SETTINGS/MCAS_config.cfg;
rm -v MCAS_v2_sim$1/SUOC/orb_in.data;
cp -v SIMSETTINGS/sim$2/MCAS_config.cfg MCAS_v2_sim$1/SETTINGS/MCAS_config.cfg;
cp -v SIMSETTINGS/sim$2/orb_in.data MCAS_v2_sim$1/SUOC/orb_in.data;
echo 'Simulation starting';
cd MCAS_v2_sim$1/;
./MCAS;
