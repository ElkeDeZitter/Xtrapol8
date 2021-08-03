import os, sys, re, numpy as np#, matplotlib.pyplot as plt
import random
"""
CLass that efficiently reads in a crystfel stream file and allow calculate statistics from the stream file as well as modifying it.
The initiat Stream.py was written by Nicolas Coquelle and later Modified by Elke De Zitter
This class should be opened with python2.7, it is not working with python3 (yet)!
"""

class Stream(object):
    def __init__(self, streamfile):
        self.streamfile = streamfile
        self.frames = []
        self.header = ''
        self.parse_stream()

    def parse_stream(self):
        append = 0
        crystal = 0
        count_shots = 0
        count_crystals = 0
        frame_stream = []
        header = 0
        indexing_methods = []

        stream = open(self.streamfile,'r')  # .readlines()

        for index, line in enumerate(stream):
            ### GEt header
            while (header == 0):
                self.header += line
                break

            ### Get beginning of an image
            if 'Begin chunk' in line:
                count_shots += 1
                append = 1
                header = 1
                # frame_stream.append(l)

            ### If
            if 'Image filename' in line: filename = line.split()[2]
            
            if 'Event: //' in line:
                event = int(re.search(r'Event:\ \/\/(.+?)$',line).group(1))
            else:
                event = 0

            if 'indexed_by' in line and 'none' not in line:
                count_crystals += 1
                crystal = 1
                frame = Frame()
                frame.indexing = line.split()[2].strip()
                if frame.indexing not in indexing_methods:
                    indexing_methods.append(frame.indexing)
                frame.filename = filename
                frame.event    = event
                try:
                    f = os.path.split(filename)[1]
                    tag = os.path.splitext(f)[0].split('tag_')[1]
                    frame.timeline = tag
                except:
                    pass
            if 'diffraction_resolution_limit' in line:
                res = float(line.split()[5])
                frame.res = res

            if 'Cell parameters' in line:
                a0, b0, c0 = line.split()[2:5]
                frame.a = float(a0)
                frame.b = float(b0)
                frame.c = float(c0)
                alpha0, beta0, gamma0 = line.split()[6:9]
                frame.alpha = float(alpha0)
                frame.beta = float(beta0)
                frame.gamma = float(gamma0)

            if "End chunk" in line:
                if crystal == 1:
                    frame_stream.append(line)
                    frame.all = frame_stream
                    self.frames.append(frame)
                append = 0
                frame_stream = []
                crystal = 0

            if append == 1: frame_stream.append(line)

            if count_shots % 1000 == 0: print '%7i frames parsed, %7i indexed frames found\r' % (count_shots, count_crystals),
            sys.stdout.flush()

        print '%7i frames parsed, %7i indexed frames found\r' % (count_shots, count_crystals),
        sys.stdout.flush()
        self.images           = count_shots
        self.indexed_images   = count_crystals
        self.indexing_methods = indexing_methods
        #self.indexed_images   = len(self.frames)
        #print 'indexing methods:',indexing_methods
        
        if 'Begin chunk' in self.header:
            self.header = re.sub('----- Begin chunk -----\n','',self.header) #quick and dirty removing the 'begin chunk' line that is added to the header
        
    
    def get_index_rate(self):
        """
        returns indexing rate = total number of images in stream / indexed images in stream
        """
        return float(self.indexed_images)/self.images
    
    def get_indexing_per_method(self):
        """
        returns the number of indexed images per indexing method
        """
        d = {}
        for method in self.indexing_methods:
            i = 0
            for f in self.frames:
                if f.indexing == method:
                    i+=1
            d[method]=i
        return d
    
    def get_total_number_of_cystals(self):
        """
        return the amount of crystals. This can be larger than the number of indexed images when multiple crystals were identified on a single image, but can never be smaller than the number of indexed images.
        """
        i = 0
        for f in self.frames:
            i+= self.get_number_of_indexed_crystals_in_frame(f.all)
        return i
    
    def get_stream_summary(self):
        """
        Print some stats about the indexing
        """
        print("number of processed images: %d" %(self.images))
        print("number of indexed images: %d" %(self.indexed_images))
        d = self.get_indexing_per_method()
        for method in d:
            print("   %s: %d" %(method, d[method]))
        print("Indexing rate: %.4f" %(self.get_index_rate()))
        print("Number of unindexed images: %d" %(self.images - self.indexed_images))
        print("Number of crystals: %d" %(self.get_total_number_of_cystals()))
    
    def get_cell_stats(self):
        """
        get statistics on cell_parameters
        modif : Paula Oeser
        """
        aas = np.array([f.a for f in self.frames])
        bbs = np.array([f.b for f in self.frames])
        ccs = np.array([f.c for f in self.frames])
        alphas = np.array([f.alpha for f in self.frames])
        betas = np.array([f.beta for f in self.frames])
        gammas = np.array([f.gamma for f in self.frames])

        aas_av = np.average(aas)
        aas_stdev = np.std(aas)
        bbs_av = np.average(bbs)
        bbs_stdev = np.std(bbs)
        ccs_av = np.average(ccs)
        ccs_stdev = np.std(ccs)
        alphas_av = np.average(alphas)
        alphas_stdev = np.std(alphas)
        betas_av = np.average(betas)
        betas_stdev = np.std(betas)
        gammas_av = np.average(gammas)
        gammas_stdev = np.std(gammas)

        return aas_av, aas, aas_stdev, bbs_av, bbs, bbs_stdev, ccs_av, ccs, ccs_stdev, alphas_av, alphas, alphas_stdev, betas_av, betas, betas_stdev, gammas_av, gammas, gammas_stdev
    
    def get_score(self):
        """
        score = indexrate/product-of-stds-on-axis
        """
        rate = self.get_index_rate()
        
        _, aas_stdev, _, bbs_stdev, _, ccs_stdev = self.get_cell_stats()
        stdev_product = aas_stdev * bbs_stdev * ccs_stdev
        
        return rate / stdev_product

    def get_image_list(self, image_list):
        """
        Split stream according to a list with image and event numbers
        """
        with open(image_list) as i:
            lst = i.read().split('\n')
        self.imagelist = []
        for fle in lst:
            self.imagelist += [f for f in self.frames if f.filename in fle and fle.endswith("//%d" %(f.event))]# and f not in self.imagelist]
            
    def get_image_lst_2(self, image_list):
        """
        Split stream according to a list with image and event numbers, but much faster than get_image_list if the list also contains images that are not present in the stream.
        """
        with open(image_list) as i:
            lst = i.read().split('\n')
        lst = np.asarray(lst)
        frame_lst = ['%s //%d'%(frame.filename, frame.event) for frame in self.frames]
        frame_lst = np.asarray(frame_lst)
        selection = np.nonzero(np.in1d(lst,frame_lst))[0]
        lst = np.take(lst, selection)
        self.imagelist = []
        for fle in lst:
            self.imagelist +=[f for f in self.frames if f.filename in fle and fle.endswith("//%d" %(f.event))]# and f not in self.imagelist]

            
    def save_image_list(self, filename):
        out = open(filename, 'w')
        print >> out, self.header,
        for f in self.imagelist:
            print >> out, ''.join(f.all),
        print 'New stream saved in file %s' % filename
        out.close()
        
    def sort(self, attribute, length=None, max_res=None):

        if attribute == 'res': return self.sortbyres(length, max_res)

        if attribute == 'a':
            if length == None:
                return sorted(self.frames, key=lambda x: x.a)
            else:
                length = min(length, self.indexed_images)
                return sorted(self.frames, key=lambda x: x.res)[0:length]

        if attribute == 'b':   return sorted(self.frames, key=lambda x: x.b)

        if attribute == 'c':   return sorted(self.frames, key=lambda x: x.c)

        if attribute == 'name': self.frmes_res = sorted(self.frames, key=lambda x: x.filename)

    def sortbyres(self, length, max_res):

        if length is None:
            if max_res is None:
                self.frames_res = sorted(self.frames,
                                         key=lambda x: x.res)  # return sorted(self.frames, key=lambda x: x.res)
            else:
                self.frames_res = sorted([i for i in self.frames if i.res <= max_res], key=lambda x: x.res)

        if max_res == None:
            length = min(length, self.indexed_images)
            self.frames_res = sorted(self.frames, key=lambda x: x.res)[0:length]
        else:
            length = min(length, len(nl))
            self.frames_res = sorted([i for i in self.frames if i.res <= max_res], key=lambda x: x.res)[:length]


    def truncate(self,start,stop,filename):
        with open(filename,'w') as f:
            f.write(self.header)
            count = start
            while count <= stop-1:
                f.write(''.join(self.frames[count].all))
                count += 1


    def save_res(filename):
        pass

    def select_indexing_methods(self, *args):

        self.indexing = []
        for meth in args:
            self.indexing += [f for f in self.frames if f.indexing == meth]

    def save_indexed(self, filename):
        if not hasattr(self, 'indexing'):
            print 'Please select the different indexing methods to be saved with "select_indexing_methods"'
            return
        else:
            out = open(filename, 'w')
            print >> out, self.header,
            for f in self.indexing:
                print >> out, ''.join(f.all),
            print 'New stream saved in file %s' % filename
        out.close()

    def save_individual(self):
        for f in self.frames:
            root = f.filename.split('.h5')[0]
            out = open(root + '.stream', 'w')
            print >> out, self.header
            print f.all
            print >> out, ''.join(f.all)

    def save_truncate(self, root, n):
        total = self.indexed_images
        chunk = total / n
        i = 0 
        num = total
        for i in range(0,n):#while i < n:
            fout = '%s_%i.stream'%(root,num)
            print 'Saving %s' %fout
            out = open(fout, 'w')
            print >> out, self.header,
            for i in range(0,num):
                print >> out, ''.join(self.frames[i].all),
            out.close()
            num -= chunk
            i += 1


    def save_filenames(self, filenames,filename):
        out = open(filename, 'w')
        print >> out, self.header,
        for frame in self.frames:
            if frame.filename in filenames: print >> out, ''.join(f.all)
        out.close()


    def select_cell(self, a, b, c, da, db, dc):
        # Add an option to only select a single axis
        # User should use 0,0 for a,da for example
#        a = a / 10.
#        b = b / 10.
#        c = c / 10.
#        da /= 10.
#        db /= 10.
#        dc /= 10.
        if not hasattr(self, 'a'): self.a = np.array([f.a for f in self.frames])
        if not hasattr(self, 'b'): self.b = np.array([f.b for f in self.frames])
        if not hasattr(self, 'c'): self.c = np.array([f.c for f in self.frames])
        print self.a
        indices_a = np.where(np.logical_and(self.a > a - da, self.a <= a + da))
        indices_b = np.where(np.logical_and(self.b[indices_a] > b - db, self.b[indices_a] <= b + db))
        indices_c = np.where(np.logical_and(self.c[indices_b] > c - dc, self.c[indices_b] <= c + dc))

        temp_indices = [indices_b[0][i] for i in indices_c[0]]
        self.final_indices = [indices_a[0][i] for i in indices_c[0]]
        print '%i crystals out of %i have been selected ' % (len(self.final_indices), self.indexed_images)
        print 'To save a new stream with these cell parameters distribution: please use save_cell("out.stream)"'

    # final_indices
    def select_cell_parameter(self, parameter, x, dx):
        """
        Select images based on a single cell parameter and spread. Cell parameters given in nm or degrees.
        """
        if parameter == 'a':
            if not hasattr(self, 'a'):
                self.x = self.a = np.array([f.a for f in self.frames])
            else:
                self.x = self.a
        elif parameter == 'b':
            if not hasattr(self, 'b'):
                self.x = self.b =np.array([f.b for f in self.frames])
            else:
                self.x = self.b
        elif parameter == 'c':
            if not hasattr(self, 'c'):
                self.x = self.c = np.array([f.c for f in self.frames])
            else:
                self.x = self.c
        elif parameter == 'alpha':
            if not hasattr(self, 'alpha'):
                self.x = self.alpha= np.array([f.alpha for f in self.frames])
            else:
                self.x = self.alpha
        elif parameter == 'beta':
            if not hasattr(self, 'beta'):
                self.x = self.beta= np.array([f.beta for f in self.frames])
            else:
                self.x = self.beta
        elif parameter == 'gamma':
            if not hasattr(self, 'gamma'):
                self.x = self.gamma= np.array([f.gamma for f in self.frames])
            else:
                self.x = self.gamma
        else:
            print("select a parameters: a, b, c, alpha, beta or gamma")
            return
        
        indices_x = np.where(np.logical_and(self.x > x - dx, self.x <= x + dx))
        self.final_indices = indices_x[0]
        print '%i crystals out of %i have been selected ' % (len(self.final_indices), self.indexed_images)
        print 'To save a new stream with these cell parameters distribution: please use save_cell("out.stream)"'


    def get_stats(self, time=False, plot=False, nbins=20, time_binning=1000):
        if time == True:
            out = open("time_stats.txt", 'w')
            threshold = time_binning
            if threshold == 1: threshold = 2
            nl = sorted(self.frames, key=lambda x: x.timeline)
            i = 0
            while ( i + threshold < len(nl)):
                t = []
                res = []
                a = []
                b = []
                c = []
                for i in range(i, i + threshold):
                    t.append(int(self.frames[i].timeline))
                    res.append(self.frames[i].res)
                    a.append(self.frames[i].a)
                    b.append(self.frames[i].b)
                    c.append(self.frames[i].c)

                t = np.array(t)
                res = np.array(res)
                a = np.array(a)
                b = np.array(b)
                c = np.array(c)
                print >> out, '%10i %7.4f %7.4f %7.4f%7.4f%7.4f %7.4f %7.4f%7.4f' % (np.average(t),
                                                                                     np.average(res),
                                                                                     np.std(res),
                                                                                     np.average(a),
                                                                                     np.std(a),
                                                                                     np.average(b),
                                                                                     np.std(b),
                                                                                     np.average(c),
                                                                                     np.std(c))

        else:
            self.res = np.array([f.res for f in self.frames])
            print np.min(self.res), np.max(self.res), np.median(self.res)
            hist_res = np.histogram(self.res, bins=nbins, range=(np.min(self.res), np.max(self.res)))
            if plot == True:
                self.plot(hist_res)
            self.a = np.array([f.a for f in self.frames])
            self.b = np.array([f.b for f in self.frames])
            self.c = np.array([f.c for f in self.frames])
            print 'cell param a - min: %6.2f - max: %6.2f - median: %6.2f' % (
                self.a.min(), self.a.max(), np.average(self.a))
            print 'cell param b - min: %6.2f - max: %6.2f - median: %6.2f' % (
                self.b.min(), self.b.max(), np.average(self.b))
            print 'cell param c - min: %6.2f - max: %6.2f - median: %6.2f' % (
                self.c.min(), self.c.max(), np.average(self.c))


    def save_cell(self, filename):
        if not hasattr(self, 'final_indices'):
            print 'Please select a new distribution of cell parameters with select_cell'
            return
        else:
            out = open(filename, 'w')
            print >> out, self.header,
            for index in self.final_indices:
                print >> out, ''.join(self.frames[index].all),
            print 'New stream saved in file %s' % filename

    def plot(self, his):
        hist, bins = his
        width = 0.7 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist, align='center', width=width)
        plt.show()

    def save_random_selection(self, n, outstreamfilename):
        total = self.indexed_images
        sele = random.sample(range(total),n)
#        fout = '%s_%iimages.stream'%(root,n)
        fout = outstreamfilename
        out = open(fout, 'w')
        print >> out, self.header, #print(self.header, file=out)
        print 'Saving %d frames to %s' %(n, fout)
        for i in sele:
            print >> out, ''.join(self.frames[i].all), #print(''.join(self.frames[i].all), file=out)
        out.close()

    def save_random_indexed_selection(self, root, n):
        if not hasattr(self, 'indexing'):
            print 'Please select the different indexing methods to be saved with "select_indexing_methods"'
            return
        else:
            total = len(self.indexing)
            sele = random.sample(range(total),n)
            fout = '%s_%iindexed.stream'%(root,n)
            out = open(fout, 'w')
            print >> out, self.header,
            print 'Saving %d indexed frames to %s' %(n, fout)
            for i,f in enumerate(self.indexing):
                if i in sele:
                    print >> out, ''.join(f.all),
                    
    def get_number_of_indexed_crystals_in_frame(self, frame_all):
        #return frame.all.count('--- Begin crystal\n') #depends too much on exact spelling
        return len([line for line in frame_all if 'Begin crystal' in line])
    
    def remove_crystal_from_multi(self, frame_all):
        begins = [i for i,line in enumerate(frame_all) if 'Begin crystal' in line]
        ends = [i for i,line in enumerate(frame_all) if 'End crystal' in line]
        del(frame_all[begins[-1]:ends[-1]+1])
        return frame_all
                    
    def save_random_indexed_crystals(self, root, n):
        if not hasattr(self, 'indexing'):
            print 'Please select the different indexing methods to be saved with "select_indexing_methods"'
            return
        else:
            total = len(self.indexing)
            sele = random.sample(range(total),n)
            len(sele)
            fout = '%s_%icrystals.stream'%(root,n)
            out = open(fout, 'w')
            print >> out, self.header,
            print 'Saving %d indexed crystals to %s' %(n, fout)
            for i,f in enumerate(self.indexing):
                if i in sele:
                    out_frame = list(f.all)
                    if self.get_number_of_indexed_crystals_in_frame(out_frame) == 0: #this should not normally not occur
                        s = next(iter(set(range(i+1, total)) - set(sele)))
                        sele.append(s)
                        print 'adding number to random selection',s
                        try:
                            assert sele.count(s) == 1
                        except AssertionError:
                            print 'random number of crystals will be smaller'
                    else:
                        while self.get_number_of_indexed_crystals_in_frame(out_frame) >1:
                            out_frame = self.remove_crystal_from_multi(out_frame)
                        else:
                            assert self.get_number_of_indexed_crystals_in_frame(out_frame) == 1, 'Number of crystal:%d for %s' %(self.get_number_of_indexed_crystals_in_frame(out_frame),f.filename)
                            print >> out, ''.join(out_frame),
           

# out=open('time_evolution.txt','w')

#for i in range(0,len(t)):
#   print >> out, "%10i %5.2f %10.7f %10.7f %10.7f"%(int(t[i]),
#   						    float(res[i]),
#						    float(a[i]),
#						    float(b[i]),
#						    float(c[i]))

class Frame():
    def init(self):
        self.a = 0
        self.b = 0
        self.c = 0
        self.alpha = 90.
        self.beta = 90.
        self.gamma = 90.
        self.res = 5.
        self.index = 'none'
        self.filename = 'example.h5'
        self.event = 0
        self.timeline = 0
        self.indexing = ''

#s=Stream('iris_nobkgsub_zaef_rings_nocen.stream')

#s=Stream('test.stream')
#nl = s.sort('res',max_res=4)

#for f in s.sort(key=lambda x: x.res)
#for f in nl: print f.res, f.filename
#s.get_stats(time=False,plot=False, time_binning=1)
#s.select_cell(8.37,9.80,14.34,0.1,0.1,0.1)
#s.save_cell('final.stream')
#s.select_indexing_methods('mosflm-comb-nolatt-nocell','dirax-comb-nolatt-nocell')

#s.save_indexed('test_indexed_stream')

