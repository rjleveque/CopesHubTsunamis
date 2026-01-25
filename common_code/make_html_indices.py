
import os,sys,glob,shutil

topdir = os.getcwd()

check_clobber = False

for (dirpath, subdirs, files) in os.walk('.'):
    files = os.listdir(os.path.abspath(dirpath))
    files.sort()
    #print('+++ dirpath = ',dirpath)
    #print('+++ subdirs = ',subdirs)
    #print('+++ files = ',files)
    os.chdir(dirpath)
    if 'index.html' in files:
        print('Directory %s already has an index' % dirpath)
	if check_clobber:
	    clobber = raw_input('Clobber it? [n] ')
	    if clobber.lower()[0] == 'y':
		shutil.move('index.html','index_backup.html')
		print('Moved to index_backup.html')

    files = os.listdir('.')
    files.sort()
    if 'index.html' not in files:
        print('Creating index.html for directory %s' % dirpath)
        with open('index.html','w') as f:
            f.write('<html>\n<h1>index of %s</h1>\n<ul>' % dirpath)
            for file in files:
                f.write('<li><a href="%s">%s</a>\n' % (file,file))
            f.write('</ul>\n</html>')
    os.chdir(topdir)
            
