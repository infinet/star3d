#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sqlite3
from math import sin, cos, pow

import matplotlib
font = {'size': 10}
matplotlib.rc('font', **font)
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

APPDIR = os.path.abspath(os.path.dirname(__file__))
DATADIR = 'data'
DB_FILE = os.path.join(APPDIR, DATADIR, 'cross.sqlite')

# IV/27A     HD-DM-GC-HR-HIP-Bayer-Flamsteed Cross Index    (Kostjuk, 2002)
# ftp://cdsarc.u-strasbg.fr/cats/IV/27A/catalog.dat
CATALOG = os.path.join(APPDIR, DATADIR, 'catalog.dat')

# I/311      Hipparcos, the New Reduction       (van Leeuwen, 2007)
# ftp://cdsarc.u-strasbg.fr/cats/I/311/hip2.dat.gz
HIP2 = os.path.join(APPDIR, DATADIR, 'hip2.dat')

# IAU constellations
CONST_CN = os.path.join(APPDIR, 'constellation-cn.txt')

PARSEC = 3.26156  # light years
COLOR_THRESH = 0  # B-V index threshold for blue and red
NAME_THRESH = 3.5  # show star name if brighter
MAX = 2000  # MAX distance of star to plot


def initdb():
    if os.path.exists(DB_FILE):
        return  # the database already generate, do nothing

    if not os.path.exists(CATALOG):
        print ("Can't find star names cross index. Please download it from "
               "ftp://cdsarc.u-strasbg.fr/cats/IV/27A/catalog.dat, "
               "save it to %s" % CATALOG)
        exit(2)

    if not os.path.exists(HIP2):
        print ("Can't find Hipparcos catalog. Please download it from "
               "ftp://cdsarc.u-strasbg.fr/cats/I/311/hip2.dat.gz, "
               "uncompress it and save to %s" % HIP2)
        exit(2)

    conn = sqlite3.connect(DB_FILE)
    sql = '''
CREATE TABLE IF NOT EXISTS hip2 (
    hip INTEGER PRIMARY KEY,
    ra TEXT,
    de TEXT,
    plx TEXT,
    hpmag TEXT,
    bv TEXT
);

CREATE TABLE IF NOT EXISTS constel (
    hd INTEGER PRIMARY KEY,
    dm TEXT,
    gc INTEGER,
    hr INTEGER,
    hip INTEGER,
    vmag REAL,
    bayer TEXT,
    cst TEXT
);

CREATE TABLE IF NOT EXISTS constel_name (
    id INTEGER PRIMARY KEY,
    fullname TEXT,
    chn TEXT,
    abr TEXT
);
'''
    conn.cursor().executescript(sql)
    conn.commit()
    read_hip2()
    read_constel()
    read_constel_fullname()


def fortran_parsefmt(fmt):
    import re
    RE_FORTRAN_FMT = re.compile(r'(\d*)([A-Z]+)(\d*)\.?\d*')
    m = RE_FORTRAN_FMT.findall(fmt)
    fw = []
    # transform Fortran format into (rep, TYPE, width)
    for x in m:
        count = int(x[0]) if x[0] else 1
        width = int(x[2]) if x[2] else 1
        fw.append((count, x[1], width))
    return fw


def fortran_read(lines, fmt):
    ''' read fortran formated file into array '''
    fw = fortran_parsefmt(fmt)
    res = []
    for line in lines:
        line.replace('D', 'E')  # replace fortran 77 double precision format
        tmp = []
        i = 0
        for rep, ftype, length in fw:
            for cnt in xrange(rep):
                field = line[i: i + length]
                i += length
                tmp.append(field)
        res.append(tmp)
    return res


def read_hip2():
    # format for Readme of I/311
    HIP2_FMT = 'I6,1X,I3,1X,I1,1X,I1,1X,F13.10,1X,F13.10,1X,F7.2,1X,F8.2,1X,F8.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,F6.2,1X,I3,1X,F5.2,1X,I2,1X,F6.1,1X,I4,1X,F7.4,1X,F6.4,1X,F5.3,1X,I1,1X,F6.3,1X,F5.3,1X,F6.3,1X,15F7.2'
    fp = open(HIP2)
    farray = fortran_read(fp, HIP2_FMT)

    conn = sqlite3.connect(DB_FILE)
    db = conn.cursor()
    sql = ('insert or replace into hip2 (hip,ra,de,plx,hpmag,bv) '
          ' VALUES (?,?,?,?,?,?)')
    # insert to db
    count = 0
    for a in farray:
        hip = int(a[0])
        ra = a[8].strip()
        de = a[10].strip()
        plx = a[12].strip()
        hpmag = a[38].strip()
        bv = a[46].strip()

        db.execute(sql, (hip, ra, de, plx, hpmag, bv))
        count += 1
    conn.commit()
    print "wrote %d records to hip2" % count


def read_constel():
    CAT_FMT = 'I6,1X,A12,1X,I5,1X,I4,1X,I6,1X,2I2,F5.2,1X,A1,2I2,F4.1,1X,F5.2,1X,I3,1X,A5,1X,A3'
    fp = open(CATALOG)
    farray = fortran_read(fp, CAT_FMT)

    conn = sqlite3.connect(DB_FILE)
    db = conn.cursor()
    sql = ('insert or replace into constel (hd,dm,gc,hr,hip,vmag,bayer,cst) '
          ' VALUES (?,?,?,?,?,?,?,?)')
    # insert to db
    count = 0
    for a in farray:
        hd = int(a[0])
        dm = a[2].strip()

        gc = int(a[4]) if a[4].strip() else None
        hr = int(a[6]) if a[6].strip() else None
        hip = int(a[8]) if a[8].strip() else None
        bayer = a[23].strip() if a[23].strip() else None

        try:
            vmag = float(a[19]) if a[19].strip() else None
        except ValueError:
            print 'valueerror'
            print a[19]
        cst = a[25].strip()
        db.execute(sql, (hd, dm, gc, hr, hip, vmag, bayer, cst))
        count += 1
    conn.commit()
    print "wrote %d records to constellation table" % count


def read_constel_fullname():
    fp = open(CONST_CN)
    conn = sqlite3.connect(DB_FILE)
    db = conn.cursor()
    sql = ('insert or replace into constel_name (abr,fullname,chn) '
          ' VALUES (?,?,?)')
    count = 0
    for line in fp:
        abr = line[:3]
        fullname = line[6: 37].strip().decode('utf8')
        chn = line[37:].strip().decode('utf8')
        db.execute(sql, (abr, fullname, chn))
        count += 1
    conn.commit()
    print "wrote %d records to constellation table" % count
    return


def sph2xyz(ra, de, r):
    z = r * sin(de)
    r_on_xy = r * cos(de)
    x = r_on_xy * cos(ra)
    y = r_on_xy * sin(ra)
    return (x, y, z)


def mag2size(mag):
    ''' calc star size based on mag'''
    s = int(pow(10, -0.3 * mag) * 100)
    return 1 if s == 0 else s


def plot_star3d(cst, mag=6.5):
    stars = find_constel(cst, mag)
    pnts = {}
    for s in stars:
        if s['mag'] > mag:
            continue
        if not s['hip']:
            continue
        if s['ly'] < 0:
            continue
        if s['ly'] > MAX:
            continue

        tmp = {}
        tmp['ly'] = s['ly']
        tmp['hip'] = s['hip']
        tmp['xyz'] = sph2xyz(s['ra'], s['de'], s['ly'])
        if s['color'] > COLOR_THRESH:
            tmp['color'] = 'r'
        else:
            tmp['color'] = 'b'
        tmp['size'] = mag2size(s['mag'])
        tmp['bayer'] = s['bayer']
        tmp['mag'] = s['mag']
        #print 'hip%d, %s, mag = %f, ly=%f, size=%d' % (s['hip'], s['bayer'],
        #                                    s['mag'], s['ly'], tmp['size'])
        pnts[s['hip']] = tmp

    print 'ploting %d objects' % len(pnts.keys())
    fig = plt.figure()
    #fig = plt.figure(figsize=(10, 7.5), dpi=75, linewidth=0.8)
    ax = fig.gca(projection='3d')
    for k, v in pnts.iteritems():
        ax.scatter(v['xyz'][0], v['xyz'][1], v['xyz'][2],
                    c=v['color'], s=v['size'])
        if v['mag'] < 3:
            ax.text(v['xyz'][0], v['xyz'][1], v['xyz'][2], v['bayer'])
        if v['ly'] > 2000:
            ax.text(v['xyz'][0], v['xyz'][1], v['xyz'][2], 'HIP%d' % v['hip'])

    # Sun
    ax.scatter(0, 0, 0, c='y', s=20)
    ax.text(0, 0, 0, 'Sun')

    ax.set_xlim(0, MAX)
    ax.set_ylim(0, MAX)
    ax.set_zlim(-1 * MAX / 2, MAX / 2)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.text2D(0.03, 0.97, '%s Stars closer than %d light years' % (cst, MAX),
            transform=ax.transAxes)
    plt.show()

    return  # do not generate output image

    STEP = 3
    for n in range(0, 170, STEP):
        angle = 190 + n
        ax.view_init(0, angle)
        plt.savefig(os.path.join(APPDIR, DATADIR, 'movie%03d.png' %
            n), bbox_inches='tight')


def find_constel(cst, mag=6.5):
    sql = ('select h.hip,h.ra,h.de,h.plx,h.hpmag,h.bv,c.bayer from '
            ' (select * from constel where cst=?) c '
            ' left outer join hip2 h '
            ' on h.hip=c.hip ')
    conn = sqlite3.connect(DB_FILE)
    db = conn.cursor()

    db.execute(sql, (cst,))
    rows = db.fetchall()
    res = []
    for row in rows:
        if not row[0]:
            continue
        tmp = {'hip': row[0], 'ra': float(row[1]), 'de': float(row[2]),
                'color': float(row[5]), 'bayer': row[6],
                'ly': PARSEC * 1000 / float(row[3]), 'mag': float(row[4])}
        if tmp['mag'] < mag and tmp['ly'] > 0:
            res.append(tmp)

    return res


def show_constel():
    conn = sqlite3.connect(DB_FILE)
    db = conn.cursor()
    sql = 'select id,fullname,chn,abr from constel_name order by id'
    db.execute(sql)
    rows = db.fetchall()

    for i in xrange(44):
        line = '%2d. %-22s %-10s  %2d. %-22s %-10s' % (
            rows[i][0], rows[i][1], rows[i][2],
            rows[i + 44][0], rows[i + 44][1], rows[i + 44][2])
        print line

    res = {}
    for r in rows:
        res[r[0]] = r[3]

    return res


def main():
    initdb()
    d = show_constel()
    sel = raw_input('Select Constellation Number: ')
    sel = int(sel)
    cst = d[sel]
    plot_star3d(cst)


if __name__ == "__main__":
    main()
