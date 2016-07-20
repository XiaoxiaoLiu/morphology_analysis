for  i in {1..100}
  do 
     python morph_cluster_on_simulation.py  -i  /data/informatics/kylel/data/morph_149n_ap/morph_149n_ap_e$i.csv  -o /data/informatics/kylel/data/morph_149n_ap/morph_149n_ap_e$i -m ap  > /data/informatics/kylel/data/morph_149n_ap/morph_149n_ap_e$i.log
     python morph_cluster_on_simulation.py  -i  /data/informatics/kylel/data/morph_149n_ward/morph_149n_ward_e$i.csv  -o /data/informatics/kylel/data/morph_149n_ward/morph_149n_ward_e$i -m ward > /data/informatics/kylel/data/morph_149n_ward/morph_149n_ward_e$i.log
done
