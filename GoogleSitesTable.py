import os
from ResolvedStellarPops import fileIO
import numpy as np
import sys
import subprocess

def get_z(filename):
    return float(os.path.split(filename)[1].split('_')[1].replace('Z',''))

def get_age(filename):
    return float(os.path.split(filename)[1].split('_')[2].replace('A','').replace('.png',''))

def trilegal_diag_table(image_location):
    """
    python GoogleSitesTable.py just run this in the directory with Marco's
    plots. 
    
    The printing at the bottom is for phil, so comment it out, change it, or
    ignore it.
    
    -- will make html for CMD and LF page with 5x5 image table (of 16 images)
    -- will make two tar.gz files with cmd*png and lf*png
    -- to change the tar commands, look at the os.system commands at line 53
    (works like perl)
    """
    image_base, header, footer, titlefmt, cellfmt, line = defaults()
    endrow = '</tr>\n'
    startrow = '<tr>'
    cell = '<td>%s</td>'
    image_garbage = '<div style="display: block; margin-right: auto; '
    image_garbage += 'margin-left: auto; text-align: center; ">'
    image_garbage += '<a imageanchor="1" href="%s%s"><img src="%s%s" border="0"'
    image_garbage += 'width="400" height="300"></a></div>'
    thismix = os.path.split(image_location)[1]
    web_image_loc = '/'.join((image_base, '%s/' % thismix))
    
    plt_types = ['cmd', 'lf']
    for plt_type in plt_types:
        outfile = os.path.join(image_location, plt_type+'.html')
        imgs = fileIO.get_files(image_location,'%s*png' % plt_type)
        zs = np.unique([get_z(z) for z in imgs])
        ages = np.unique([get_age(a) for a in imgs])

        content = "".join(['<th>Z = %g </th>' % z for z in zs])
        line = '%s<th width=60px>&nbsp;</th>%s%s' % (startrow, content, endrow)

        for age in ages:
            line += '%s<th>%g Gyr</th>' % (startrow, age)
            for z in zs:
                img = '_'.join((plt_type,'Z%.2e' % z, 'A%.2e.png' % age))
                line += cell % (image_garbage % (web_image_loc, img,
                                                 web_image_loc, img))
            line += endrow
        
        quick_write(outfile, header, line, footer)
        
        here = os.getcwd()
        os.chdir(image_location)

        os.system('tar -cvf %ss.tar %s*png' % (plt_type, plt_type))
        os.system('gzip %ss.tar' % plt_type)
        os.chdir(here)
        
    if 'philrose' in image_location:
        print 'scp %s*.gz philrose@portal.astro.washington.edu:/www/astro/users/philrose/html/Research/tpagbcalib/.'%fmt
        print 'ssh philrose@portal.astro.washington.edu'
        print 'cd /www/astro/users/philrose/html/Research/tpagbcalib'
        print 'mkdir %s'%thismix
        print 'mv *gz %s'%thismix

def defaults():
    '''
    gotta be a better way...
    these are just globals in a way. Most are googlesites computer vomit.
    '''
    image_base = 'http://www.astro.washington.edu/users/philrose/Research/tpagbcalib'
    header = '<table border="1" bordercolor="#888888" style="border-color:rgb(136,136,136);border-collapse:collapse"><tbody>\n'
    footer = '</tbody></table>\n'
    titlefmt = '<tr><td style=\'text-align: center;\'><h3>%s</h3></td></tr>\n'
    cellfmt = '<tr><td><div style=\'display:block;text-align:center;margin-right:auto\'><a href=\'%s\'><img height=\'%i\' src=\'%s\' width=\'%i\'></a><br></div></td></tr>\n'
    line = ''
    return image_base, header, footer, titlefmt, cellfmt, line
    
def w_h_frompng(img):
    '''
    finds the width had height of a png.
    '''
    info = subprocess.Popen([r"file",img],stdout=subprocess.PIPE).communicate()[0]
    wh = info.split(',')[1]
    w,h = wh.split('x')
    w = int(w.strip())
    h = int(h.strip())
    return w,h
     
def side_by_side(local_dir,image_dir,outfile,
                 extra_left='diag.',extra_right='diag_',zip_it=False):
    '''
    two column table images set by width and height of left image. Nrows is 
    the number of pngs in the directory. 
    input
    local_dir: str where images are living
    image_dir: directory on the server (http etc is given in defaults)
    outfile: html to write to
    extra_left: file locator *extra_left*png for objects in column 1
    extra_right: same as left but for column 2
    zip_it: make a tar.gz ball.    
    
    '''
    image_base,header,footer,titlefmt,cellfmt,line = defaults()
    
    abs_image_dir = '/'.join((image_base,os.path.join(image_dir,'diag')))
    imgsL = fileIO.get_files(local_dir,'*%s*png'%extra_left)
    imgsR = fileIO.get_files(local_dir,'*%s*png'%extra_right)
    
    titlefmt = titlefmt.replace('<td','<td colspan="2"')   
    newrow = '<tr>'
    endrow = '</tr>\n'
    cellfmt = cellfmt.replace(newrow,'').replace(endrow,'')
    
    for imgL,imgR in zip(imgsL,imgsR):
        imgtitle = os.path.split(imgL)[1].split('_')[0]
        img_locL = '/'.join((abs_image_dir,os.path.split(imgL)[1]))
        img_locR = '/'.join((abs_image_dir,os.path.split(imgR)[1]))
        [fileIO.ensure_file(img) for img in (imgL,imgR)]
        w,h = w_h_frompng(imgL)
        line += titlefmt % imgtitle
        line += newrow
        line += cellfmt % (img_locL,h/2,img_locL,w/2)
        line += cellfmt % (img_locR,h/2,img_locR,w/2)
        line += endrow
        
    quick_write(outfile,header,line,footer)
    if zip_it==True: zip_em(local_dir,image_dir,extra='*diag')
    return

def zip_em(local_dir,image_dir,extra=''):
    '''
    create a tar ball of the images to put on server.
    '''
    os.system('tar -cvf %s.tar %s/%s*png'%(local_dir,local_dir,extra))
    os.system('gzip %s.tar'%local_dir)
    if 'philrose' in image_dir:
        print 'scp *.gz philrose@portal.astro.washington.edu:/www/astro/users/philrose/html/Research/tpagbcalib/.'
        print 'ssh philrose@portal.astro.washington.edu'
        print 'cd /www/astro/users/philrose/html/Research/tpagbcalib'
        print 'mkdir %s'%os.path.split(image_dir)[1]
        print 'mv *gz %s'%os.path.split(image_dir)[1]
    return
    
def quick_write(outfile,header,line,footer):
    '''
    write header, line, and footer to a file
    '''
    out = open(outfile,'w')
    out.write(header)
    out.write(line)
    out.write(footer)
    out.close()
    print 'wrote %s'%outfile
    return
    
def one_col(image_location, local_dir='.', outfile='pngs.list', extra='', zip_it=True):
    '''
    one column separated by titles
    input
    local_dir: str where images are living
    image_dir: directory on the server (http etc is given in defaults)
    outfile: html to write to
    extra: file locator *extra*png
    zip_it: make a tar.gz ball.    
    '''
    
    image_base, header, footer, titlefmt, cellfmt, line = defaults()
    
    abs_image_dir = '/'.join((image_base, image_location))
    imgs = fileIO.get_files(local_dir, '*%s*png' % extra)
    
    for img in imgs:
        imgtitle = os.path.split(img)[1].split('_')[0]
        img_loc = '/'.join((abs_image_dir, os.path.split(img)[1]))
        w,h = w_h_frompng(img)
        line += titlefmt % imgtitle
        line += cellfmt % (img_loc, h, img_loc, w)
    
    quick_write(outfile, header, line, footer)
    if zip_it is True:
        zip_em(local_dir, image_location)
    return
    

if __name__ == "__main__":
    image_location = sys.argv[1:]
    one_col(*image_location)
    
        
