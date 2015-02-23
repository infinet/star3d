#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import ftplib
import gzip
import hashlib
import sqlite3
from math import sin, cos, pow

APPDIR = os.path.abspath(os.path.dirname(__file__))
DB_DIR = os.path.join(APPDIR, 'data')
DB_FILE = os.path.join(APPDIR, 'star3d.sqlite')

# IV/27A     HD-DM-GC-HR-HIP-Bayer-Flamsteed Cross Index    (Kostjuk, 2002)
# I/311      Hipparcos, the New Reduction       (van Leeuwen, 2007)
urls = {
'crossindex': {'url': 'ftp://cdsarc.u-strasbg.fr/cats/IV/27A/catalog.dat',
               'desc': 'HD-DM-GC-HR-HIP-Bayer-Flamsteed Cross Index',
               'sha1sum': 'e7771854cacf410ffb3065811295755f1aaffd70',
               'local': os.path.join(DB_DIR, 'catalog.dat')},
'hip2':  {'url': 'ftp://cdsarc.u-strasbg.fr/cats/I/311/hip2.dat.gz',
          'desc': 'Hipparcos',
          'sha1sum': '854e7a65f40f7e7beb8c54f3ba60de9ab2a8cb54',
          'local': os.path.join(DB_DIR, 'hip2.dat.gz')}
}

# IAU constellations
CONST_CN = os.path.join(APPDIR, 'constellation-cn.txt')

PARSEC = 3.26156  # light years
COLOR_THRESH = 0  # B-V index threshold for blue and red
NAME_THRESH = 3.5  # show star name if brighter


def ftp_get(url, output):
    if url.startswith('ftp://'):
        tmp = url[6:].split('/')
    else:
        tmp = url.split('/')
    host = tmp[0]
    fname = tmp[-1]
    path = '/'.join(tmp[1:-1])
    ftp = ftplib.FTP(host)
    ftp.login()
    ftp.cwd(path)
    with open(output, 'wb') as fp:
        ftp.retrbinary('RETR %s' % fname, fp.write)
    ftp.quit()


def verify_file(fname, sha1):
    ''' verify file's sha1sum, return True when match'''
    h = hashlib.new('sha1')
    content = open(fname).read()
    h.update(content)
    if sha1 == h.hexdigest():
        return True
    else:
        return False


def initdb():
    if os.path.exists(DB_FILE):
        return  # the database already generate, do nothing

    if not os.path.exists(DB_DIR):
        os.mkdir(DB_DIR)

    for k,v in urls.iteritems():
        if not os.path.exists(v['local']):
            print "Can't find %s. Begin download from %s" % (v['desc'],
                                                             v['url'])
            ftp_get(v['url'], v['local'])

        if not verify_file(v['local'], v['sha1sum']):
            print 'sha1sum does NOT match, remove %s and try again' % v['local']
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
    target = urls['hip2']
    print 'Reading %s into database' % target['desc']
    fp = gzip.open(target['local'])
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
    target = urls['crossindex']
    print 'Reading %s into database' % target['desc']
    fp = open(target['local'])
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


def plot_star3d(cst, mag=6.5, output='screen', max_dist=1000):
    '''
    Args:
        cst: constellation IAU abbreviation
        mag: magnitude
        output: choose between screen and image
        max_dist: maximum distance of star to plot, default 1000 light years
    '''

    import matplotlib
    font = {'size': 10}
    matplotlib.rc('font', **font)
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    if output == 'screen':
        fig = plt.figure()
    else:
        matplotlib.use('Agg')
        fig = plt.figure(figsize=(10, 7.5), dpi=75, linewidth=0.8)

    stars = find_constel(cst, mag)
    pnts = {}
    for s in stars:
        if s['mag'] > mag:
            continue
        if not s['hip']:
            continue
        if s['ly'] < 0:
            continue
        if s['ly'] > max_dist:
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

    print 'Plotting %d objects' % len(pnts.keys())
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

    ax.set_xlim(0, max_dist)
    ax.set_ylim(0, max_dist)
    ax.set_zlim(-1 * max_dist / 2, max_dist / 2)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.text2D(0.03, 0.97, '%s Stars closer than %d light years' % (cst, max_dist),
            transform=ax.transAxes)

    if output == 'screen':
        plt.show()
        return  # do not generate output image

    STEP = 10
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


def show_constel(chn = False):
    conn = sqlite3.connect(DB_FILE)
    db = conn.cursor()
    sql = 'select id,fullname,chn,abr from constel_name order by id'
    db.execute(sql)
    rows = db.fetchall()

    if chn:
        for i in xrange(44):
            line = '%2d. %-22s %-10s  %2d. %-22s %-10s' % (
                rows[i][0], rows[i][1], rows[i][2],
                rows[i + 44][0], rows[i + 44][1], rows[i + 44][2])
            print line
    else:
        tmp = []
        for i in xrange(88):
            if (i + 1) % 3 == 0:
                tmp.append('%2d. %-22s' % (rows[i][0], rows[i][1]))
                print ' '.join(tmp)
                tmp = []
            else:
                tmp.append('%2d. %-22s' % (rows[i][0], rows[i][1]))

        print ' '.join(tmp)

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

    sel = raw_input('Limit maximum distance of stars to [Default '
            '1000 light years]: ')
    try:
        dist = int(sel)
    except ValueError:
        dist = 1000
    print 'Plotting stars closer than %d light years' % dist
    plot_star3d(cst, mag=6.5, output='screen', max_dist=dist)


if __name__ == "__main__":
    main()
