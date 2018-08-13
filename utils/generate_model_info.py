
import os
import glob


CM_LIBRARY = '/rna/auto-traveler/data/crw-cm'


all_cm = 'all.cm'  # file with all CMs
all_cm_path = os.path.join(CM_LIBRARY, all_cm)

if not os.path.exists(all_cm_path):
    cmd = 'cat {CM_LIBRARY}/*.cm > {CM_LIBRARY}/all.cm'.format(CM_LIBRARY=CM_LIBRARY)
    os.system(cmd)

with open(os.path.join(CM_LIBRARY, 'modelinfo.txt'), 'w') as f:
    line = '*all*    -    -    %s\n' % all_cm
    f.write(line)
    for cm in glob.glob('%s/*.cm' % CM_LIBRARY):
        if all_cm in cm:
            continue
        model_name = os.path.basename(cm).replace('.cm', '')
        line = "%s    SSU    Bacteria    %s\n" % (model_name, os.path.basename(cm))
        f.write(line)
print 'Done'
