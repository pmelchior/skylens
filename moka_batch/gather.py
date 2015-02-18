#!/bin/env python

import sqlite3
import glob
import os, sys
import time

if len(sys.argv) < 2:
    print "usage: " + argv[0] + " <output DB>"
    exit(0)

outfile = sys.argv[1]
conn = None
cur = None
if os.path.exists(outfile):
    conn = sqlite3.connect(outfile) # master db: must exists with table
    cur = conn.cursor()

dirs = glob.glob('moka_batch_*')
now = time.time()
for dir in dirs:
    dbfile = dir + "/SHEAR_ACC.db"
    # check if file exists and is older than 10 minutes
    # avoids opening an active db and break the transaction lock
    if os.path.exists(dbfile) and now - os.path.getmtime(dbfile) > 10*60:
        print dbfile
        
        # check if output DB is present
        # need proper table for INSERT to be successful
        if conn is not None:
            connX = sqlite3.connect(dbfile)
            curX = connX.cursor()
        
            # copy contents of dbfile to master
            cur.execute("ATTACH '" + dbfile + "' as dbX")
            cur.execute("INSERT INTO shear_accuracy SELECT * from dbX.shear_accuracy")
            cur.execute("DETACH DATABASE dbX")
        # if not: move first suitable dbfile to outfile and open connection
        else:
            os.system("mv " + dbfile + " " + outfile)
            conn = sqlite3.connect(outfile) # master db: must exists with table
            cur = conn.cursor()
        
        # clean up
        os.system('rm -rf ' + dir)



