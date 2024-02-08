import os
import sqlite3
from sqlite3 import OperationalError
import glob

og_dbs = glob.glob("../input/gtex_v8_sqtl_dbs_mashr/*.db")
num_og = len(og_dbs)

con = sqlite3.connect("../input/combine_sqtl.db")
cur = con.cursor()
print("Connected to main db")

num_done = 0
for index, single_db in enumerate(og_dbs):
    db_fname = os.path.basename(single_db)
    db_tiss_name = db_fname.split('mashr_')[1].split('.')[0] # Depends on mashr_<tissue_name>.db format
    
    cur.execute("ATTACH '{}' as db{}".format(single_db, index))
    # print("Attached {}".format(single_db))

    combine = "INSERT INTO " + "extra" + " SELECT * FROM db{}.".format(index) + "extra"
    # print(combine)
    cur.execute(combine)

    combine = "INSERT INTO " + "weights" + " SELECT '{tiss}' ||'.'|| gene, rsid, '{tiss}' ||'.'|| varID, ref_allele, eff_allele, weight FROM db{ind}.".format(tiss=db_tiss_name, ind=index) + "weights"
    print(combine)
    cur.execute(combine)

    con.commit()

    detach = False

    cur.execute("detach database db{}".format(index))

cur.close()
con.close()
