#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.

"""Helper functions for testing error correcting code functions"""
import time, gc
from codinglib.util import print_progress, weight, flatten_once
from codinglib.code import DecodingFailedError

# A timing function which keeps a global state
_timeping = None
def time_ping(mes = None, done=False):
    """Convenient time keeping during execution. Call time_ping() to reset
    timer, call time_ping(mes="...") to print a sub-time since last call to
    time_ping naming the time span with mes. Call time_ping(done=True) to print a
    final message summing up the total time spent since time_ping() was
    called."""
    global _timeping
    newtime = time.time()
    if mes or done:
        if mes:
            print "%s took %.2f s" % (mes, newtime-_timeping[1])
        if done:
            print "Total time was %.2f s" % (newtime-_timeping[0])
            _timeping = None
        else:
            _timeping[1] = newtime
    else:
        _timeping = [newtime,newtime]


# Structures for holding the timing results
__clock_call_run = 0
__clock_call_run_d = dict()
__clock_call_dict = { __clock_call_run: __clock_call_run_d }
__clock_call_stopped = False
__clock_call_overview_title = "clock call overview"
def clock_call(f):
    """Decorator for collecting timing statistics on functions. Timings are
    collected simply as elapsed time between start and end of all functions. For
    a given function, we remember such an elapsed duration for every call to
    that function. Furthermore, "runs" are supported which allows for logically
    separating such monitoring between logically independent runs of the code
    (for instance, in order to calculate a mean expenditure).
    
    Due to the cost of collecting these statistics, the timings are skewed from
    reality by a factor which increases with the frequency with which function
    is called. Monitored "outer" functions calling many sub-functions which are
    being monitored will be less reliable.

    The statistics can be manipulated and accessed through clock_call_new_run,
    clock_call_reset, clock_call_results, clock_call_details.
    """

    def timed(*args, **kw):
        if __clock_call_stopped:
            return f(*args, **kw)
        else:
            before = time.time()
            result = f(*args, **kw)

            if not f.__name__ in __clock_call_run_d:
                __clock_call_run_d[f.__name__] = []
            __clock_call_run_d[f.__name__].append(time.time() - before)
            return result

    return timed

def clock_call_new_run():
    """Begin a new "run" of the code for which the collected timings will be
    separated from already collected timings."""
    global __clock_call_run_d
    global __clock_call_stopped
    global __clock_call_run
    if __clock_call_run_d:
        if not __clock_call_stopped:
            __clock_call_run_d["end"] = time.time()
        __clock_call_run += 1
        __clock_call_run_d = dict()
        __clock_call_dict[__clock_call_run] = __clock_call_run_d
        __clock_call_stopped = False
    __clock_call_run_d["start"] = time.time()

def clock_call_stop_run():
    """Stop the current "run" of the code. It is only necessary to call
    clock_call_stop_run before clock_call_new_run if one is interested in
    omitting operations and calls from being monitored between the two runs."""
    global __clock_call_stopped
    if not __clock_call_stopped:
        __clock_call_run_d["end"] = time.time()
        __clock_call_stopped = True

def clock_call_reset():
    """Reset all collected clock_call data"""
    global __clock_call_run
    global __clock_call_run_d
    global __clock_call_stopped
    global __clock_call_dict
    __clock_call_run = 0
    __clock_call_run_d = dict()
    __clock_call_dict = { __clock_call_run: __clock_call_run_d }
    __clock_call_stopped = False

def clock_call_results():
    """Return a sensible summary of the collected timings"""
    all_funcs = set()
    for run in __clock_call_dict:
        all_funcs = all_funcs.union(list(__clock_call_dict[run]))
    all_funcs.remove("start")
    all_funcs.remove("end")
    results = dict()
    def add(d, key, add):
        if key in d:
            d[key] += add
        else:
            d[key] = copy(add)
    overview = dict()
    results[__clock_call_overview_title] = overview
    for run in __clock_call_dict:
        run_d = __clock_call_dict[run]
        add(overview, "total_runs", 1)
        if "start" in run_d:
            run_time = run_d["end"] - run_d["start"]
            add(overview, "total_run_time", run_time)
            add(overview, "all_run_times", [run_time])
        clocked_time = sum(sum(run_d[f]) if f in run_d else 0 for f in all_funcs) # Time caught in clock_calls
        add(overview, "clocked_time", clocked_time)

    if "total_run_time" in overview:
        overview["avg_run_time"] = overview["total_run_time"]/overview["total_runs"]
        overview["median_run_time"] = median(overview["all_run_times"])
        overview["fraction_clocked_time"] = overview["clocked_time"]/overview["total_run_time"]
    for f in all_funcs:
        f_d = dict()
        results[f] = f_d
        for run in __clock_call_dict:
            run_d = __clock_call_dict[run]
            if not f in run_d:
                add(f_d, "absent", 1)
            else:
                calls = run_d[f]
                time_in_run = sum(calls)
                add(f_d, "present_in_runs", 1)
                add(f_d, "total_calls", len(calls))
                add(f_d, "total_time", time_in_run)
                add(f_d, "all_calls", calls)
                add(f_d, "all_run_times", [time_in_run])
                if "start" in run_d:
                    add(f_d, "all_run_percentages", [time_in_run/(run_d["end"] - run_d["start"])])
        f_d["avg_time_in_run"] = f_d["total_time"]/f_d["present_in_runs"]
        f_d["median_time_in_run"] = median(f_d["all_run_times"])
        f_d["avg_time_in_call"] = f_d["total_time"]/f_d["total_calls"]
        f_d["median_time_in_call"] = median(f_d["all_calls"])
        f_d["rel_std_dev_time_in_call"] = 100*std(f_d["all_calls"])/f_d["avg_time_in_call"] if f_d["avg_time_in_call"] > .001 and len(f_d["all_calls"]) > 1 else -1
        if "all_run_percentages" in f_d:
            f_d["median_run_percentage"] = 100*median(f_d["all_run_percentages"])
        f_d.pop("all_calls")
        f_d.pop("all_run_times")
        f_d.pop("all_run_percentages")
    return results

def clock_call_print_results():
    """Print a sensible summary of the collected timings"""
    print "\nTIMING RESULTS\n"
    results = clock_call_results()
    name_width = int(30)
    count = ("%6i ", int(7))
    secs  = ("% 9.3f secs.", int(15))
    perc = ("% 8.3f%% ", int(10))
    overview = results[__clock_call_overview_title]
    if overview["total_runs"] > 1:
        print(("No. of runs:\t" + count[0]) % overview["total_runs"])
    if "total_run_time" in overview:
        print(("Total run time:\t" + secs[0]) % overview["total_run_time"])
        print(("Avg. run time:\t" + secs[0]) % overview["avg_run_time"])
        print(("Med. run time:\t" + secs[0]) % overview["median_run_time"])
    columns = [ ("total_calls", "#calls", count),
                ("total_time", "Total time", secs),
                ("median_time_in_run", "Med. time/run", secs),
                ("median_time_in_call", "Med. time/call", secs),
                ("rel_std_dev_time_in_call", "stddev/avg", perc),
                ("median_run_percentage", "% of run", perc) ]
    pstr = "%*s" + "|%*s"*len(columns)
    # Print table headline
    print(pstr % tuple([name_width, "Function/Method"] + list(flatten_once( (typ[1], cname) for (_,cname,typ) in columns ))))
    print("-"*(name_width + sum( 1+typ[1] for (_,_,typ) in columns)))
    # Print functions list
    for fun in results:
        if fun != __clock_call_overview_title:
            fun_d = results[fun]
            pstr = "%*s|" + "|".join(cstr[0] for (_,_,cstr) in columns)
            print(pstr % tuple([name_width, fun]+[ fun_d[key] if key in fun_d else "N/A" for (key,_,_) in columns ]))
    
    tab_to_last_col = " "*(name_width + 1 + sum(columns[i][2][1]+1 for i in range(len(columns)-1)))
    # Print % of run summation
    print(tab_to_last_col + "-"*columns[-1][2][1])
    print(tab_to_last_col + (perc[0] % (100*overview["fraction_clocked_time"])))

    

def clock_call_details(fun_name):
    """Return all collected timings for the named function"""
    details = dict()
    for run in __clock_call_dict:
        if fun_name in __clock_call_dict[run]:
            details[run] = __clock_call_dict[run][fun_name]
    return details
            
    
    

####################################################################
#          RANDOM CONSTRUCTIONS
####################################################################

def random_error_pos(n, errs):
    """Return precisely errs different random numbers between 0 and n-1. errs
must be <= n"""
    return sample(range(n), errs)
    
def random_error_vec_at_pos(n, F, errPos):
    """Construct a random error vector with the given error positions (other
positions are zero)"""
    vec = [F.zero()]*n
    for i in errPos:
        while vec[i].is_zero():
            vec[i] = F.random_element()
    return vector(vec)

def random_error(n, F, errs):
    """Construct a random error vector over the field F of length n with
precisely errs non-zero entries."""
    return random_error_vec_at_pos(n, F, random_error_pos(n, errs))

def error_pos(err_vec):
    """Return the non-zero positions of the given error vector"""
    n = len(err_vec)
    return [ i for i in range(0, n) if not err_vec[i].is_zero() ]

####################################################################
#          TEST DECODERS
####################################################################

def decoding_instance(C, nerrs):
    """Return a decoding instance (r,e,c) for hard-decision decoding exactly nerrs errors"""
    c = C.random_codeword()
    e = random_error(C.n, C.F, nerrs)
    r = c + e
    return (r, e, c)

class DecodingTestFailed(Exception):
    """Exception to throw when a decoding test failed.
    Includes a partial report until the failure."""
    def __init__(self, mes, report):
        self.mes = mes
        self.report = report
    def __str__(self):
        return self.mes

def test_decoder(C, decoder, N=1, nerrs=None, nerrs_max=None, nerrs_min=1, continue_on_error=False, silent=False, garbage_collect=False, progress=False):
    r"""Test a decoding algorithm on a linear block code with a set of N random
    codewords, subjected to random errors patterns, supposedly decodable.

    Returns a report being a list of trial results. Each trial result is a tuple
    of the form `success * time [* list size ]`, where the last element is
    always 1 for non-list decoders. The test is immediately halted on the first
    decoding failure if `continue_on_error=False`.
    More information is also printed out, except if silent=True.

    For list decoders, decoding success is defined as the correct codeword being on the returned list.

    Takes numerous optional arguments:

    - N                 - Number of received words to decode. Default = 1.
    - nerrs_max         - Maximal number of errors to add to the codeword. Default = decoder.decoding_radius().
    - nerrs_min         - Minimal number of errors to add to the codeword. Default = 1.
    - nerrs             - If set, denotes the exact number of errors to add to the codeword. If set, then
                          nerrs_max and nerrs_min are ignored.
    - silent            - Do not print out any information.
    - continue_on_error - On a decoding failure or error, return immediately if
                          continue_on_error is False. Otherwise, continue the test
    - garbage_collect   - Force garbage collection between each decoded word.
    - progress          - Report progress during the test
    """
    import gc
    report = []
    def add_to_report(nerrs, success, elapsed, ans):
        if decoder.is_list_decoder():
            report.append((nerrs, success, elapsed, len(ans)))
        else:
            report.append((nerrs, success, elapsed, 1))
    def myprint(s):
        if not silent:
            print(s)
    if nerrs:
        nerrs_max = nerrs
        nerrs_min = nerrs
    if nerrs_max == None:
        nerrs_max = decoder.decoding_radius()
    progress_size = N//100 if N >= 100 else 1
    for iters in range(N):
        nerrs = randint(nerrs_min, nerrs_max)
        (r,e,c) = decoding_instance(C, nerrs)
        before = time.time()
        try:
            ans = decoder.decode(r)
            fail = None
        except DecodingFailedError as exc:
            fail = "DECODING FAILURE"
            ans = []
        elapsed = time.time() - before
        if garbage_collect:
            gc.collect()
        if fail is None and c!=ans and (decoder.is_list_decoder() and not (c in ans)):
            fail = "ERRONEOUS DECODING"
        if fail is None:
            add_to_report(nerrs, True, elapsed, ans)
        else:
            add_to_report(nerrs, False, elapsed, ans)
            if not continue_on_error:
                raise DecodingTestFailed("%s IN %sth trial (decoder: %s)\nr = %s!" % (fail, iters+1, decoder, r), report)
        if progress and iters % progress_size == 0:
            print "%s%% done" % (round(iters/N*100)+1)
    times = [ elapsed for _,_,elapsed,_ in report ]
    times.sort()
    successes = sum(1 for _,succ,_,_ in report if succ)
    err_dict = dict()
    for nerrs,_,_,_ in report:
        if nerrs in err_dict:
            err_dict[nerrs] += 1
        else:
            err_dict[nerrs] = 1
    myprint("%s out of %s decoding(s) successful. no. tests for each no. errs: %s\nDecoding took on median %7.3f secs. (from %7.3f to %7.3f secs.)"\
            % (successes, N, err_dict, times[N//2], times[0], times[-1]))
    if decoder.is_list_decoder():
        list_sizes = dict()
        for _,_,_,ell in report:
            if ell in list_sizes:
                list_sizes[ell] += 1
            else:
                list_sizes[ell] = 1
        myprint("Encountered output list sizes: %s" % list_sizes)
    return report

def test_decoded_information(f, C, r, tau):
    """Returns whether the information word f, when turned into a codeword,
    is within distance tau from the received word r"""
    return weight(vector(r) - vector(C.encode(f))) <= tau
