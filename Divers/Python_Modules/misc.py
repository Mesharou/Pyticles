"""
Miscellaneous things that are too small to warrant a module.
"""

from string import Template

import time, os, errno, re

# This module must not call logutils.getLogger upon import.
# The problem is that Bunch is imported by system.__init__.
# If getLogger were called here, then a root logger would
# be initialized before any other module can make what is
# supposed to be the *first* call to getLogger.


# based on Robert Kern's Bunch:
class Bunch(dict):
    """
    A dictionary that also provides access via attributes.

    Additional methods update_values and update_None provide
    control over whether new keys are added to the dictionary
    when updating, and whether an attempt to add a new key is
    ignored or raises a KeyError.

    The Bunch also prints slightly differently than a normal
    dictionary, using str() instead of repr() for its values,
    and in key-sorted order.
    """
    def __init__(self, *args, **kwargs):
        """
        *args* can be dictionaries, bunches, or sequences of
        key,value tuples.  *kwargs* can be used to initialize
        or add key, value pairs.
        """
        dict.__init__(self)
        self.__dict__ = self
        for arg in args:
            self.update(arg)
        self.update(kwargs)

    def __str__(self):
        items = list(self.items())
        items.sort()
        s = "{%s}" % (' '.join(["%s: %s" % (repr(k),str(v)) for (k,v) in items]))
        ## fix the formatting later
        #slist = ['Dictionary with access to the following as attributes:']
        #slist.append('\n'.join(self.keys()))
        #return '\n'.join(slist) + '\n'
        return s

    def update_values(self, *args, **kw):
        """
        arguments are dictionary-like; if present, they act as
        additional sources of kwargs, with the actual kwargs
        taking precedence.

        One reserved optional kwarg is "strict".  If present and
        True, then any attempt to update with keys that are not
        already in the Bunch instance will raise a KeyError.
        """
        strict = kw.pop("strict", False)
        newkw = dict()
        for d in args:
             newkw.update(d)
        newkw.update(kw)
        self._check_strict(strict, newkw)
        dsub = dict([(k, v) for (k, v) in list(newkw.items()) if k in self])
        self.update(dsub)

    def update_None(self, *args, **kw):
        """
        Similar to update_values, except that an existing value
        will be updated only if it is None.
        """
        strict = kw.pop("strict", False)
        newkw = dict()
        for d in args:
             newkw.update(d)
        newkw.update(kw)
        self._check_strict(strict, newkw)
        dsub = dict([(k, v) for (k, v) in list(newkw.items())
                                if k in self and self[k] is None])
        self.update(dsub)

    def _check_strict(self, strict, kw):
        if strict:
            bad = set(kw.keys()) - set(self.keys())
            if bad:
                bk = list(bad)
                bk.sort()
                ek = list(self.keys())
                ek.sort()
                raise KeyError(
                    "Update keys %s don't match existing keys %s" % (bk, ek))


class Cachefile(object):
    """
    class for cache file manilulation: init, read, write
    file is key, value; anything after a '#' is ignored
    """

    def __init__(self, cachefile=None, contents='parameters'):
        if cachefile is None:
            raise IOError('no file specified')

        self.cachefile = cachefile
        self.contents = contents
        self.comments = []

    #-------
    def init(self, *args, **kw):
        '''
        start the file, eg metafile with information to cache,
                           or which is not in the database
        initialize file with *args (dict, bunch, tuples) and **kw
        '''
        # only init if empty or nonexistent
        if os.path.exists(self.cachefile):
            from pycurrents.system import logutils
            # This late and local initialization is needed; otherwise
            # a root logger is created when logutils itself is
            # imported, which defeats its functionality.
            log = logutils.getLogger(__file__)
            log.warning('file %s exists.  not initializing' %
                          (self.cachefile))
        self.cachedict=Bunch({})
        self.cachedict.update(*args, **kw)

        # initialize the comments here
        self.comments = ['#this file was automatically generated. DO NOT EDIT',
               '#\n'
               '# written %s\n' % (nowstr()),
               '# this file contains ' + self.contents + '\n',
               '#name, value pairs',
               '#--------------------',
              ]

        self.write()
    #-------
    def _strip_comments(self, lines):
        '''
        strip comments, return a dictionary
        comment protocol: put "#" at the front of line or whitespace split chunk
        '''
        cache_dict = {}
        comments = []
        for line in lines:
            line.rstrip()
            parts = line.split()
            for ipp in range(len(parts)):
                if parts[ipp][0] == '#':   #simple comment strip
                    comments.append(' '.join(parts[ipp:]))
                    parts=parts[:ipp]
                    break
            if len(parts) == 2:
                cache_dict[parts[0]] = parts[1]


        return cache_dict, comments

    #-------
    def add_comments(self, lines):
        '''
        add comments
        '''
        comments = []
        for line in self.comments:
            if line not in comments:
                comments.append(line)

        for line in lines:
            if line not in self.comments:
                self.comments.append(line)

    #------
    @staticmethod
    def _parse(arg):
        try:
            ret = float(arg)
            if int(ret) == ret:
                return int(ret)
            return ret
        except ValueError:
            if arg.lower() == 'none':
                return None
            return arg

    #----
    def _translate(self):
        '''
        translate  'None' into None, read numbers as float or int
        '''
        for key, val in list(self.cachedict.items()):
            self.cachedict[key] = self._parse(val)

    #-------
    def read(self):
        '''
        metafile with information to cache, or which is not in the database
        read; sets attribute "cache" (Bunch, from file contents)

        '''
        ## This needs help ("intnames" is clunky)

        # read if present
        if os.path.exists(self.cachefile) is False or (
                             os.path.exists(self.cachefile) is True and
                             os.path.getsize(self.cachefile)==0):
            raise IOError('cannot read cachefile %s', self.cachefile)

        lines = open(self.cachefile,'r').readlines()
        cdict, self.comments = self._strip_comments(lines)
        self.cachedict = Bunch(cdict)
        self._translate()

    #-------
    def write(self, *args, **kw):
        '''
        update variables in self.cachedict from *args and  kw.keys
        key must exist in cache *before* calling write method
        then write to cachefile

        '''

        # read anyway, initialize it if necessary

        self.cachedict.update_values(*args, **kw)

        llist=[]
        kk= list(self.cachedict.keys())
        kk.sort()
        for k in kk:
            llist.append('%20s   %s' %(k,self.cachedict[k]))
        cf = open(self.cachefile,'w')
        cf.write('\n'.join(self.comments))
        cf.write('\n')
        cf.write('\n'.join(llist))
        cf.write('\n')
        cf.close()


def sleep(secs):
    '''
    This works just like time.sleep except that it traps an error.
    The IOError "unknown error 514" seems to be caused by a linux
    kernel bug.
    '''
    try:
        time.sleep(secs)
    except IOError:
        pass


def guess_comment(fname, cstr='%#', numlines=10):
    '''
    guess the comment character from a line of text
    (matlab processing uses '%', python processing uses '#')

    returns first comment character match
    '''
    F=open(fname,'r')
    comment = None
    for num in range(numlines):
        line = F.readline().strip()
        for cc in cstr:
            if len(line) > 0:
                if cc == line.split()[0][0]:
                    return cc


def nowstr():
    '''
    get utc time from computer clock, return string "yyyy/mm/dd hh:mm:ss"
    '''
    return time.strftime("%Y/%m/%d %H:%M:%S")


class ScripterBase(object):
    """
    Base class for classes that write standalone scripts.

    The script is left behind for later command-line use.
    The code in the script is run using exec.

    This is designed for quick_adcp.

    """

    script_head = ""
    script_body = ""
    script_tail = ""
    defaultparams = {}

    def __init__(self, *args, **kw):
        """
        The optional first argument is a dictionary of options.
        Additional options may be specified via keyword arguments.
        """
        if len(args):
            self.optdict = args[0]
        else:
            self.optdict = {}
        self.kw = kw
        fname = self.__class__.__name__ + '_script.py'
        self.script_filename = kw.get('script_filename', fname)
        self.params = self.defaultparams.copy()
        for k, v in list(self.optdict.items()):
            if k in self.params:
                self.params[k] = v
        self.params.update(kw)

    def __call__(self, **kw):
        """
        Call can be used for stand-alone plotting, with any kw
        params supplied here overriding those from instantiation.

        Subclass should call this method at the start of its own
        __call__ method.
        """
        self.params.update(kw)
        self.process_params()
        # Turn all parameters into attributes, for convenience.
        for k, v in self.params.items():
            setattr(self, k, v)
        # Actual plotting code follows this in subclass.

    def process_params(self):
        """
        Subclass will use this for validation etc.
        """
        pass

    def fill_strings(self):
        self.head = Template(self.script_head).substitute(self.params)
        self.body = Template(self.script_body).substitute(self.params)
        self.tail = Template(self.script_tail).substitute(self.params)

    def write(self, filename=None):
        self.process_params()
        self.fill_strings()
        if filename is None:
            filename = self.script_filename
        script = '\n'.join([self.head, self.body, self.tail, ''])
        open(filename, 'w').write(script)

    def run(self):
        self.process_params()
        self.fill_strings()
        exec(self.body)

def safe_makedirs(*args, **kw):
    """
    Replacement for os.makedirs.

    os.makedirs raises an exception if the requested path
    already exists; safe_makedirs raises
    an exception only if the requested path cannot be made, or
    if it exists but is not a directory.
    """
    try:
        os.makedirs(*args, **kw)
    except OSError as e:
        if e.errno != errno.EEXIST or not os.path.isdir(args[0]):
            raise

