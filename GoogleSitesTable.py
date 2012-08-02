# GoogleSitesTable.py
"""
python GoogleSitesTable.py just run this in the directory with Marco's plots. 

The printing at the bottom is for phil, so comment it out, change it, or ignore it.

-- will make html for CMD and LF page with 5x5 image table (of 16 images)
-- will make two tar.gz files with cmd*png and lf*png
-- to change the tar commands, look at the os.system commands at line 53 (works like perl)
"""
import os
import glob
import numpy as np
import sys

def get_z(filename):
    return float(filename.split('_')[1].replace('Z',''))

def get_age(filename):
    return float(filename.split('_')[2].replace('A','').replace('.png',''))

def main(image_location):
    header = '<table style="border-collapse: collapse; border-top-color: rgb(136, 136, 136); border-right-color: rgb(136, 136, 136); border-bottom-color: rgb(136, 136, 136); border-left-color: rgb(136, 136, 136); border-top-width: 1px; border-right-width: 1px; border-bottom-width: 1px; border-left-width: 1px; " border="1" bordercolor="#888" cellspacing="0">\n<tbody>\n'
    footer = '</tbody></table>\n'
    endrow = '</tr>\n'
    startrow = '<tr>'
    cell = '<td>%s</td>'
    image_garbage = '<div style="display: block; margin-right: auto; margin-left: auto; text-align: center; "><a imageanchor="1" href="%s%s"><img src="%s%s" border="0" width="400" height="300"></a></div>'
    thismix = os.path.split(os.getcwd())[1]
    image_loc = os.path.join(image_location,'%s/'%thismix)
    
    plt_types = ['cmd','lf']
    for plt_type in plt_types:
        outfile = '_'.join((os.path.split(os.getcwd())[1],plt_type+'.html'))
        out = open(outfile,'w')
        imgs = glob.glob1('.','%s*png'%plt_type)
        zs = np.unique([get_z(z) for z in imgs])
        ages = np.unique([get_age(a) for a in imgs])
        
        line = startrow+'<th width=60px>&nbsp;</th>'+"".join(['<th>Z = %g </th>'%z for z in zs])+endrow
        for age in ages:
            line+= startrow+'<th>%g Gyr</th>'%age
            for z in zs:
                img = '_'.join((plt_type,'Z%.2e'%z,'A%.2e.png'%age))
                line+= cell%(image_garbage%(image_loc,img,image_loc,img))
            line += endrow
        
        out.write(header)
        out.write(line)
        out.write(footer)
        
        os.system('tar -cvf %ss.tar %s*png'%(plt_type,plt_type))
        os.system('gzip %ss.tar'%plt_type)
    out.close()
    if 'philrose' in image_location:
        print 'scp *.gz philrose@portal.astro.washington.edu:/www/astro/users/philrose/html/Research/tpagbcalib/.'
        print 'ssh philrose@portal.astro.washington.edu'
        print 'cd /www/astro/users/philrose/html/Research/tpagbcalib'
        print 'mkdir %s'%thismix
        print 'mv *gz %s'%thismix

if __name__ == "__main__":
    image_location = sys.argv[1]
    main(image_location)    
    
        
