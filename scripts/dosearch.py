#!/usr/bin/python

import os, subprocess, sys, getopt, glob, time, re
from time import sleep
from multiprocessing import Process

# A few useful functions
def command_output(command):
    try:
        # Note that subprocess.check_output retuns a byte string
        pipe = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        out, err = pipe.communicate()
        return out
    except:
        print('Error in {0}\n'.format(command)) 
        sys.exit(1)

def add_holes(holes, treeholes, directory, boxfile):
    command = '{0} -r {1} \'{2}\''.format(treeholes, directory, boxfile)
    byte_string = command_output(command)
    if not byte_string: return
    byte_string = byte_string.replace('root','')
    if len(byte_string) == 0 : byte_string = 'root'
    new_holes = set(byte_string.rstrip().split('\n'))
    holes |= set(h for h in new_holes if len(h) < 121)

def add_holes_from_file(holes,fp) :
    try:
        for line in open(fp):
            hole = line.rstrip()
            if hole[0] == '1' or\
               hole[0] == '0' :
                if len(hole) < 81 :  
                    holes.add(hole)
    except:
        print('Error loading holes file {0}\n'.format(fp))

def add_words(words, fp):
    try:
        for line in open(fp):
            word = line.rstrip()
            if word[0].isdigit() or\
               word[0] == 'X'    or\
               word[0] == 'H' : continue
            else:
                if '(' in word :
                    word = re.findall('\((.*?)\)', word)[-1]
                    words.add(word)
    except:
        print('Error loading words file {0}\n'.format(fp))
        # sys.exit(1)

def run_refine(command, dest_dir) :
    pid = os.getpid()
    pid_file = dest_dir + '/' + str(pid) + '.pid'
    with open(pid_file,'a') as fp : 
        fp.write(command + '\n')
        fp.close()
    returnCode = subprocess.call(command, shell=True)
    if returnCode == 0:
        with open(pid_file,'a') as fp :
            fp.write('\ncompleted')
        sys.exit(0)
    else:
        with open(pid_file,'a') as fp :
            fp.write('\nfailed')
        sys.exit(1)

if __name__ == '__main__' :

    try:
        opts, args = getopt.getopt(sys.argv[1:],'w:p:c:d:h:r:i:t:s:',['words=','powers=','child_limit=','depth_limit=','holes=','refine=','invent_depth=','truncate_depth','word_search_depth='])
    except getopt.GetoptError as err:
        print str(err)
        print('Usage: dosearch [-w,--words <words_file>] [-p,--powers <powers_file>] [-c,--child_limit <limit>] [-d,--depth_limit <limit>] [-h,-holes <holes_file>] src_dir dest_dir')
        sys.exit(2)

    if len(args) != 2:
        print('Usage: dosearch [-w,--words <wordsfile>] [-p,--powers <powersFile>] [-c,--child_limit <limit>] [-d,--depth_limit <limit>] [-h,-holes <holesfile>] src_dir dest_dir')
        sys.exit(2)

    # Executables
    treecat = './treecat'
    treeholes = './treecat --open_holes'
    treecheck = './treecat --mark -s'
    refine = './refine_marg'

    # Set up the rest of the arguments
    src_dir = args[0]
    dest_dir = args[1]
    childLimit = 8
    depth_limit = 330

    max_size = '5000000'
    max_depth = '330'
    truncate_depth = '6'
    invent_depth = '42'
    word_search_depth = '6'
    fill_holes = ''
    improve_tree = ''
    powers_file = 'none'
    words_file = '/home/ayarmola/margulis_search/words'

    # Get config
    holes_file = None
    seen_words = set()
    for opt, val in opts:
        if opt in ('-w', '--words'):
            words_file = val
        if opt in ('-p', '--powers'):
            powers_file = val
        if opt in ('-c', '--child_limit'):
            childLimit = int(val)
        if opt in ('-d', '--depth_limit'):
            depth_limit = int(val)
        if opt in ('-h', '--holes'):
            holes_file = val
        if opt in ('-r', '--refine'):
	    refine = val
        if opt in ('-i', '--invent_depth'):
	    invent_depth = str(int(val))
        if opt in ('-t', '--truncate_depth'):
	    truncate_depth = str(int(val))
        if opt in ('-s', '--word_search_depth'):
	    word_search_depth = str(int(val))
            
    add_words(seen_words, words_file)

    # Check for incomplete trees
    subprocess.call('{0} -r {1} \'{2}\''.format(treecheck, dest_dir, ''), shell=True)

    # Get holes. Note, treecat will check that all files are complete trees
    holes = set();
    if holes_file :
        add_holes_from_file(holes, holes_file)
    else :
        add_holes(holes, treeholes, dest_dir, '')

    # Get done words
    done = set()
    try:
        done = set([os.path.basename(boxfile).replace('.out','') for boxfile in glob.glob(dest_dir + '/*.out')])
    except:
        print('Error reading {0}\n'.format(dest_dir)) 
        sys.exit(1)

    print "Launching Refine"

    # Launch the refine runs
    active_pid_to_hole = {};
    refine_run_count = 0
    child_count = 0
    wait_for_holes = False
    failed_holes = set()
    while True:
        sleep(0.01) # We don't need to to run the main loop to death since we aren't using os.wait
        open_holes = holes - done
        if len(open_holes) == 0 and refine_run_count == 0 and len(done) == 0:
            best_hole = 'root'
        else : 
            best_hole = '1'*400
        for hole in open_holes:
            if len(hole) < len(best_hole) :
                best_hole = hole    

        if len(best_hole) > depth_limit:
            if child_count > 0 :
                wait_for_holes = True
            else :
                # We only break if we don't have any more refine processes running
                break
        else :
            wait_for_holes = False

        # We now check for completed refine processes.
        if child_count >= childLimit or (child_count > 0 and len(open_holes) == 0) or wait_for_holes:
            iter_dict = dict(active_pid_to_hole)
            for done_pid, done_hole in iter_dict.iteritems() :
                pid_file = dest_dir + '/' + str(done_pid) + '.pid'
                status = command_output('tail -1 {0}'.format(pid_file))
                if 'completed' in status :
                    # We should check the output either way to make sure it is clean 
                    subprocess.call('{0} {1} \'{2}\''.format(treecheck, dest_dir, done_hole), shell=True)

                    print 'Completed {0} {1}\n'.format(done_hole,done_pid)
                    add_holes(holes, treeholes, dest_dir, done_hole)

                    num_patched = command_output('grep -c Patched {0}/{1}.err; exit 0'.format(dest_dir, done_hole)).rstrip()
                    num_unpatched = command_output('grep -c Unpatched {0}/{1}.err; exit 0'.format(dest_dir, done_hole)).rstrip()
                    hum_holes = command_output('grep -c HOLE {0}/{1}.err; exit 0'.format(dest_dir, done_hole)).rstrip()
                    
                    print 'Holes: {0} patched, {1} unpatched, {2} open holes\n'.format(num_patched, num_unpatched, int(hum_holes))

                    box_words = set()
                    add_words(box_words, '{0}/{1}.out'.format(dest_dir, done_hole))        
                    new_words = box_words - seen_words
                    seen_words |= new_words

                    bad_holes = command_output('grep HOLE {0}/{1}.err | cut -d " " -f 2; exit 0'.format(dest_dir, done_hole)).rstrip().split('\n')
                    failed_holes.update(bad_holes)

                    if len(new_words) > 0: 
                        f = open(words_file, 'a')
                        for word in new_words:
                            print 'Adding word {0}'.format(word)
                            f.write('(' + word + ')' + '\n')
                        f.close()

                    child_count -= 1
                    del active_pid_to_hole[done_pid]
                    os.remove(pid_file)
                    continue

                elif 'failed' in status :
                    # We should check the output either way to make sure it is clean 
                    subprocess.call('{0} {1} \'{2}\''.format(treecheck, dest_dir, done_hole), shell=True)
                    # If there was an error refining
                    print 'Error with pid {0}\n'.format(done_pid)
                    print 'Error refining hole {0}\n'.format(done_hole)
                    done.remove(done_hole)
                    child_count -= 1
                    del active_pid_to_hole[done_pid]
                    os.remove(pid_file)
                    continue
                else :
                    continue
            sleep(0.01) # We don't need to to run the main loop to death since we aren't using os.wait
            continue        

        # If we make it here. We are running refine
        print 'Open hole count: {0}\n'.format(len(open_holes))
        print 'Best hole: {0}\n'.format(best_hole)
        if len(failed_holes) > 0:
          print 'Deepest failed hole: {}\n'.format(sorted(failed_holes, key=len)[-1])
        else:
          print 'Deepest failed hole: None\n'

        out = dest_dir + '/' + best_hole + '.out'
        err = dest_dir + '/' + best_hole + '.err'

        if best_hole == 'root':
            pid_word_search_depth = '-1'
        else: 
            pid_word_search_depth = word_search_depth

        treecat_command = '{0} {1} {2}'.format(treecat, src_dir, best_hole)
        refine_command = refine + \
                    fill_holes + \
                    improve_tree + \
                    ' --box ' + best_hole + \
                    ' --max_depth ' + max_depth + \
                    ' --truncate_depth ' + truncate_depth + \
                    ' --invent_depth ' + invent_depth + \
                    ' --max_size ' + max_size + \
                    ' --words ' + words_file + \
                    ' --word_search_depth ' + pid_word_search_depth + \
                    ' --powers ' + powers_file + \
                    ' > ' + out  + ' 2> ' + err

        first_command = treecat_command + ' | head -1'
        first = command_output(first_command).rstrip()

        if first[:1] == 'H': # HOLE
            treecat_command = 'echo 1'

        command = treecat_command + ' | ' + refine_command
        print 'Running with run count {1}: {0}\n'.format(command, refine_run_count)
        refine_run = Process(target=run_refine, args=(command, dest_dir,))
        refine_run.start()
        pid = refine_run.pid

        child_count += 1   
        refine_run_count += 1
        done.add(best_hole)
        active_pid_to_hole[pid] = best_hole
        done_failed = set()
        for h in failed_holes:
          if h.startswith(best_hole):
            done_failed.add(h)
        failed_holes.difference_update(done_failed)
