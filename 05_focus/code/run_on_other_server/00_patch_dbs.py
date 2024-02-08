# This code replaces the rsid column of the db with the varID
# This helps the db when imported to FOCUS format since our LD
# data uses varID ~ chr<num>_<pos>_<a1>_<a2> to identify variants.
import os
import sqlite3
from sqlite3 import OperationalError
import glob

og_dbs = glob.glob("*.db")
num_og = len(og_dbs)

# con = sqlite3.connect("combine_sqtl.db")
# cur = con.cursor()
# print("Connected to main db")

num_done = 0
for index, single_db in enumerate(og_dbs):
    os.system(f"cp template_db/empty_template.db patched/{single_db}")
    db_fname = os.path.basename(single_db)
    db_tiss_name = db_fname.split('mashr_')[1].split('.')[0] # Depends on mashr_<tissue_name>.db format

    con = sqlite3.connect(f"patched/{single_db}")
    cur = con.cursor()
    
    cur.execute("ATTACH '{}' as db{};".format(single_db, index))
    # print("Attached {}".format(single_db))

    combine = "INSERT INTO " + "extra" + " SELECT * FROM db{}.".format(index) + "extra"
    # print(combine)
    cur.execute(combine)

    combine = "INSERT INTO " + "weights" + " SELECT gene, varID, varID, ref_allele, eff_allele, weight FROM db{ind}.".format(ind=index) + "weights;"
    print(db_tiss_name)
    print(combine)
    print()
    cur.executescript(combine)

    cur.execute("detach database db{};".format(index))
    cur.close()

    con.close()
