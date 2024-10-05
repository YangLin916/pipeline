from turtle import end_fill
from matplotlib.pyplot import close
from pipeline.python.pipeline import experiment 
import datajoint as dj
acq = dj.create_virtual_module('acq','acq')
import h5py
import os 
import numpy as np
from pipeline.python.pipeline.utils import h5

from commons.python.commons import lab
from stimuli.python.stimulus.stimulus import Sync,Trial


schema = dj.schema('pipeline_stimulus', create_tables=False)

@schema 
class EphysSync(dj.Computed):
    definition = """ # syncing scanimage frame times to stimulus clock
    -> experiment.Scan
    ---
    beh2stim_slope              : double
    beh2stim_intercept          : double
    npixel2beh_intercept        : double
    npixel2beh_slope            : double
    num_npixel_samples          : int 
    sync_ts=CURRENT_TIMESTAMP   : timestamp
    """

    @property 
    def key_source(self):
        return experiment.Scan() & experiment.ScanEphysLink()

    def make(self,key):
        beh2stim_intercept,beh2stim_slope,beh_fs_correction_factor = behavior_to_stimulus(key)
        num_npixel_samples,npixel2beh_intercept,npixel2beh_slope = npixel_to_behavior(key, beh_fs_correction_factor)
        
        self.insert1({**key, 'beh2stim_slope':beh2stim_slope
                      ,'beh2stim_intercept':beh2stim_intercept
                      ,'npixel2beh_intercept':npixel2beh_intercept
                      ,'npixel2beh_slope':npixel2beh_slope
                      ,'num_npixel_samples': num_npixel_samples})
        


def npixel_samples_in_stimulus_clock(npixel_samples,beh2stim_slope,beh2stim_intercept,npixel2beh_slope,npixel2beh_intercept): 
    np_in_stim_clock = beh2stim_slope*(npixel2beh_slope*npixel_samples + npixel2beh_intercept) + beh2stim_intercept
    return np_in_stim_clock



def get_real_times(relvar):
    # Convert hardware counter to real times accounting for wraparound
    #   accurate_time = get_real_times(relvar) converts counter values to real times
    #   in ms for the tuples in relvar self.

    # Assuming `relvar` is a list of tuples or a structured array

    [count, timestamper_time, session_start_time] = relvar.fetch('count', 'timestamper_time', 'session_start_time')

    timestamper_time = np.array(timestamper_time, dtype=np.float64)
    session_start_time = np.array(session_start_time, dtype=np.float64)

    # Rescale to times
    counter_rate = 10e6 / 1000  # pulses / ms (should be stored somewhere)
    counter_period = 2**32 / counter_rate

    count_time = count / counter_rate
    approximate_session_time = timestamper_time - session_start_time

    # Compute expected counter value based on CPU time
    approximate_session_periods = np.floor(approximate_session_time / counter_period)
    approximate_residual_period = np.mod(approximate_session_time, counter_period)

    # Correct edge cases where number of periods is off by one
    idx = np.where((approximate_residual_period - count_time) > counter_period / 2)
    approximate_session_periods[idx] += 1

    accurate_time = count_time + approximate_session_periods * counter_period

    return accurate_time
    
    
def get_positive_edge_time_stamps(key):
    """
    Retrieve positive edge timestamps for a given key.

    return: List of timestamps (in ms) for positive edges
    """
    try:
        # Fetch ephys keys based on the key
        ephys_keys = acq.Ephys & (experiment.ScanEphysLink & key)
        if not ephys_keys:
            raise ValueError('No ephys recording available for this session')

        # Pick the last ephys key
        if len(ephys_keys) != 1:
            raise ValueError('Multiple ephys sessions are found in one session manager session, unable to analyze')

        # Fetch ephys stop time
        ephys_start_time, ephys_stop_time = ephys_keys.fetch('ephys_start_time', 'ephys_stop_time')

        # Fetch the channel for NPSync
        channel_key = acq.TimestampSources & 'source= "NPSync"' & ephys_keys

        # Form the key expression to get timestamps
        relvar = acq.SessionTimestamps & ephys_keys & channel_key & f'timestamper_time>{ephys_start_time[0]} AND timestamper_time<{ephys_stop_time[0]}'

        # Get the times of positive edges after accounting for wraparounds
        ts = get_real_times(relvar)

    except Exception as e:
        print(f'Sync Timestamps are Unavailable: {e}')
        ts = []

    return ts
 

def get_neuropixel(scan_key):
    '''
    
    
    '''
    head_directory = '/Volumes/scratch11/Ephys/' ### this needs to be manually changed --> needs to beautomatic from the database
    experiment= dj.create_virtual_module('experiment','pipeline_experiment')
    filename = (acq.Ephys() & (experiment.ScanEphysLink() & scan_key)).fetch1('ephys_path').split('/')[2:]
    proper_filename = "/".join(filename)
    scan_filename = head_directory + "/" + proper_filename
    print(scan_filename)
    scan = h5py.File(scan_filename,'r',driver='family', memb_size=2147483648)
    sync_signal_from_ephys_file = scan['data_2'][384,:]
    Fs = scan.attrs['Fs'][0]
    scan.close
    print(f'length of neuropixel signal: {len(sync_signal_from_ephys_file)}')

    return {'npixel_scan':sync_signal_from_ephys_file, 'npixel_Fs':Fs}


def get_behavior_file(scan_key):
    '''
    

    '''
    beh_path = (experiment.Session() & scan_key).fetch1('behavior_path')
    beh_path = lab.Paths().get_local_path(beh_path)
    filename = (experiment.Scan().BehaviorFile() & scan_key).fetch1('filename')
    
    behavior_path = (experiment.Session() & scan_key).fetch1('behavior_path')
    local_path = lab.Paths().get_local_path(behavior_path)
    filename = (experiment.Scan.BehaviorFile() & scan_key).fetch1('filename')
    full_filename = os.path.join(local_path, filename)

    behfile = h5py.File(full_filename,'r',driver='family', memb_size=2147483648)
    beh_fs = behfile.attrs['AS_Fs'][0] #behavior sampling rate
    behfile.close

    data = h5.read_behavior_file(full_filename)
    

    return data, beh_fs



def behavior_to_stimulus(scan_key):

    data, beh_fs = get_behavior_file(scan_key)  
    print('inside beh to stim' , len(data), beh_fs)
    print(f'length of scanImage: {len(data["scanImage"])}')

    # Get counter timestamps and convert to seconds
    photodiode_times = h5.ts2sec(data['ts'], is_packeted=True)
    print(f'print photodiode times length: {len(photodiode_times)}')

    # Detect rising edges in scanimage clock signal (start of each frame)
    binarized_signal = data['scanImage'] > 2.7 # TTL voltage low/high threshold
    rising_edges = np.where(np.diff(binarized_signal.astype(int)) > 0)[0]

    # Get positive edges from database clock
    positive_edge_time_stamp = get_positive_edge_time_stamps(scan_key)

    # check that length of edges is the same betweeen behavior and database
    if rising_edges.shape != positive_edge_time_stamp.shape:
        print(f'edges mismatch: {len(rising_edges) - len(positive_edge_time_stamp)}')

    # Calculate behavior duration from the rising edges of the behavior file
    #beh_fs = 10000 ## you should get it from the behavior h5 files
    beh_duration = (rising_edges[-1]-rising_edges[0])/beh_fs

    # Calculate the behavior from the database in msec
    db_duration = (positive_edge_time_stamp[-1]-positive_edge_time_stamp[0])/1000

    # to correct the photodiode times 
    fs_correction_factor = beh_duration/db_duration
    photodiode_times = photodiode_times/fs_correction_factor
        
    # the intercept and slope will transform behavior times to stimulus time
    # slope is expected to be close to 1
    intercept, slope, sfts, tfts, high_residual = Sync._sync(data['syncPd'], photodiode_times, Trial() & scan_key)
    
    return intercept, slope, fs_correction_factor
    


def change_range(x,y):

    y_range = (y.max() - y.min())
    x_range = (x.max() - x.min())
    shifted_y = (((y - y.min())*x_range)/y_range)+x.min()
    return shifted_y


def get_rising_edges(x,cutoff):
    binarized = x > cutoff
    rising_edges = np.where(np.diff(binarized.astype(int)) == 1)[0]
    return rising_edges

  
# we are obtaing a conversion factor that converts neuropixel sample time referenced from 0 and using neuropixels sampling rate
# to behavior time which is time recorded by a counter
# for now, we are not using the high resolution time stamps and given that sync variability in mouse pipeline is so high, 
# maybe this high resolution timing may not be needed
def npixel_to_behavior(scan_key, beh_fs_correction_factor):
    print('inside npx to beh')
    neuropixel = get_neuropixel(scan_key)
    print(len(neuropixel))
    behavior, beh_fs = get_behavior_file(scan_key)
    print('read h5 files', beh_fs)
    npixel_sync = neuropixel['npixel_scan']
    behavior_sync = behavior['scanImage']
    npixel_Fs = neuropixel['npixel_Fs']

    print('convert beh tmst')

    behavior_sample_time = h5.ts2sec(behavior['ts'], is_packeted=True) 
    # time correction--> divide by factor
    behavior_sample_time = behavior_sample_time/beh_fs_correction_factor

    # rising edge sample indices in sync signal recorded in behavior file
    behavior_rising_edges = get_rising_edges(behavior_sync,3)

    # convert rising edge sample indices from behavior file's sync signal to time in secs
    behavior_rising_edges_ts = behavior_sample_time[behavior_rising_edges]
    baseBehTime = behavior_rising_edges_ts[0]
    behavior_rising_edges_ts_rel_zero = behavior_rising_edges_ts - baseBehTime

    #shifted_neuropixel = change_range(behavior_sync,npixel_sync)
    npixel_rising = get_rising_edges(npixel_sync,0.7)

    npx_duration = (npixel_rising[-1]-npixel_rising[0])/npixel_Fs

    # Get positive edges from database clock
    positive_edge_time_stamp = get_positive_edge_time_stamps(scan_key)

    # check that length of edges is the same betweeen behavior and database
    if len(npixel_rising) != len(positive_edge_time_stamp):
        print(f'edges mismatch: {len(npixel_rising) - len(positive_edge_time_stamp)}')
        ### make them equal

    # Calculate the duration from the database in msec
    db_duration = (positive_edge_time_stamp[-1]-positive_edge_time_stamp[0])/1000

    npx_fs_cor_factor = npx_duration/db_duration

    # rising edge sample indices in sync signal recorded in neuropixels ephys file to time in secs rel zero being first sample
    npixel_rising_ts = npixel_rising / (npixel_Fs*npx_fs_cor_factor)
    baseNpixelTime = npixel_rising_ts[0]
    npixel_rising_ts = npixel_rising_ts - baseNpixelTime

    # computed correlations between the edge times from the two files (one shifted w.r.t to another) 
    # and find the shift that yields maximum correlation
    shift_min = -10 # edges, may need to check this range
    shift_max = 10 # edges
    mismatch = [0 for i in range(shift_max-shift_min+1)]
    shift_array = [0 for i in range(shift_max-shift_min+1)]

    iter = 0
    for shift in range(shift_min, shift_max):

        # shift the appropriate vector elements
        if shift > 0:
            v1 = npixel_rising_ts[shift:]
            v2 = behavior_rising_edges_ts_rel_zero
        elif shift < 0:
            v1 = npixel_rising_ts
            v2 = behavior_rising_edges_ts_rel_zero[shift:]
        else:
            v1 = npixel_rising_ts
            v2 = behavior_rising_edges_ts_rel_zero

        # match the sizes of the vectors
        if len(v1) > len(v2):
            v1 = v1[:len(v2)]
        elif len(v1) < len(v2):
            v2 = v2[:len(v1)]

        # just compute the sum of difference squared as a measure of time mismatch in the two vectors
        mismatch[iter] = np.sum((v1 - v2) ** 2)
        shift_array[iter] = shift
        print(mismatch[iter])
        print(shift)
        iter = iter + 1

    # determine the shift for which the mismatch was minimum
    min_mismatch_idx = np.where(mismatch == np.amin(mismatch))
    if len(min_mismatch_idx) != 1:
        print("There are more than one minima in the mismatch function")   
    min_shift = shift_array[min_mismatch_idx[0][0]]
    print(f'min_shift : {min_shift}')

    # extract the proper section from the two vectors to do a regression, sync signal edges from ephys file is in samples, 
    # and from behavior file will be in seconds
    if min_shift > 0:
        v1 = npixel_rising_ts[min_shift:] # for ephys file's sync signal, we want to use sample index for regressor
        v2 = behavior_rising_edges_ts
    elif min_shift < 0:
        v1 = npixel_rising_ts
        v2 = behavior_rising_edges_ts[min_shift:]
    else:
        v1 = npixel_rising_ts
        v2 = behavior_rising_edges_ts
       
    # match the sizes of the vectors
    if len(v1) > len(v2):
        v1 = v1[:len(v2)]
    elif len(v1) < len(v2):
        v2 = v2[:len(v1)]
   
    ## linear model fitting 
    from sklearn.linear_model import TheilSenRegressor
    npixel2behavior = TheilSenRegressor() 
    npixel2behavior.fit(np.array(v1).reshape(-1,1),np.array(v2).reshape(-1,1))
    return len(npixel_sync), npixel2behavior.intercept_, npixel2behavior.coef_[0] 

    