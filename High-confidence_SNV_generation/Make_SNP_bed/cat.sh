#!/bin/bash

for i in {X,Y,MT}
do
  cat Compara.H.sap_P.tro_lastz_net_on_H.sap.$i.*.coords > Compara.H.sap_P.tro_lastz_net_on_H.sap.$i.ALL.maf.coords
  echo "Finished processing chromosome $i."
done

for i in {1..22}
do
  cat Compara.H.sap_P.tro_lastz_net_on_H.sap.$i.*.coords > Compara.H.sap_P.tro_lastz_net_on_H.sap.$i.ALL.maf.coords
  echo "Finished processing chromosome $i."
done



