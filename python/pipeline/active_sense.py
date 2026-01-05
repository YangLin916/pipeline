import datajoint as dj

from pipeline import experiment, mice


schema = dj.schema("active_sense", locals(), create_tables=True)


@schema
class MouseInfo(dj.Manual):
    definition = """ # baseline identity and restriction metadata for acrive_sense mice

    -> mice.Mice
    ---
    sex='U'                               : enum('M', 'F', 'U')        # M/F/Unknown for quick reference on the cage card
    dob = null                            : date                      # date of birth
    name = null                           : varchar(64)               # optional nickname/alias
    line=''                               : varchar(128)              # strain or line name
    genotype=''                           : varchar(255)              # genotype details
    single_housing_date = null            : date                      # start of single housing (baseline days 1-3)
    restriction_start_date = null         : date                      # start of water restriction (days 4+)
    baseline_weight_g = null              : decimal(6,3)              # average baseline weight (g)
    baseline_weight_date = null           : date                      # date baseline weight was measured
    notes=''                              : varchar(1024)             # free-text notes
    """


@schema
class DailyLog(dj.Manual):
    definition = """ # daily weight and water intake log for water-restricted mice

    -> MouseInfo
    log_date                            : date                          # date of measurement
    ---
    phase='single_housing'              : enum('single_housing', 'water_restriction', 'post_restriction')
    -> experiment.Person                # experimenter performing measurements
    body_weight_g                       : decimal(6,3)                  # mouse weight (g)
    weight_pct_baseline = null          : decimal(5,2)                  # percent of baseline body weight
    weight_delta_from_baseline_g = null : decimal(6,3)                  # weight change vs baseline (current - baseline)
    bottle_weight_before_g = null       : decimal(7,3)                  # Drinko bottle weight before adding water (retrieval weight)
    bottle_weight_after_g = null        : decimal(7,3)                  # Drinko bottle weight after adding water (refill weight)
    water_added_ml = null               : decimal(6,3)                  # amount of water added (ml) to bottle
    
    # Task specific
    task = null                         : varchar(128)                  # task name (e.g. SoundDetection) or null if no task
    task_duration_min = null            : smallint                      # duration of task in minutes
    total_trials = null                 : smallint                      # total trials performed
    correct_trials = null               : smallint                      # number of correct trials
    correct_rate = null                 : decimal(5,2)                  # correct / total (0-100 or 0-1)
    reward_size_ul = null               : decimal(5,2)                  # reward size in microliters
    task_water_ml = 0                   : decimal(6,3)                  # water consumed during task
    
    remain_water_ml = null              : decimal(6,3)                  # target_water - task_water (how much more needed)
    
    water_consumed_ml = null            : decimal(6,3)                  # calculated intake (bottle delta + task_water)
    target_water_ml = null              : decimal(6,3)                  # allowance based on baseline weight * 50 ml/kg
    minimum_water_ml = null             : decimal(6,3)                  # minimum allowance based on 25 ml/kg rule
    removed_from_restriction = 0        : boolean                       # mark if mouse taken off restriction today
    health_status = '5'                 : enum('1', '2', '3', '4', '5') # health score (1=Critical, 5=Normal)
    notes = ''                          : varchar(1024)                 # free-text notes (spillage, handling, etc.)
    """

    @staticmethod
    def compute_daily_log(animal_id, log_date, body_weight_g,
                          bottle_weight_before_g, bottle_weight_after_g,
                          water_added_ml, 
                          task=None, task_duration_min=None, 
                          total_trials=None, correct_trials=None, 
                          reward_size_ul=None, task_water_ml=0,
                          username='', health_status='5', notes=''):
        """
        Helper method to compute derived fields for a daily log entry.
        Returns a dictionary suitable for insert1.
        """
        # Fetch baseline info
        baseline = (MouseInfo & {'animal_id': animal_id}).fetch1()
        baseline_weight = float(baseline['baseline_weight_g'])

        # Calculate weight stats
        weight_pct_baseline = (body_weight_g / baseline_weight) * 100
        weight_delta_from_baseline_g = body_weight_g - baseline_weight

        # Calculate water allowances
        target_water_ml = baseline_weight * 0.05
        minimum_water_ml = baseline_weight * 0.025

        # Task Calculations
        correct_rate = None
        remain_water_ml = None
        
        # Ensure task_water_ml is float for calculation
        if task_water_ml is None:
            task_water_ml = 0.0
        else:
            task_water_ml = float(task_water_ml)
            
        if task and total_trials and total_trials > 0:
            if correct_trials is not None:
                correct_rate = float(correct_trials) / float(total_trials)
        
        # Remain water = Target - Task Water
        remain_water_ml = target_water_ml - task_water_ml

        # Water Consumed can only be calculated the NEXT day (retrospective)
        # So for the current log, it remains None.
        water_consumed_ml = None
        
        # Handle None conversions for DB
        if correct_rate is None: correct_rate = 0.0
        
        return dict(
            animal_id=animal_id,
            log_date=log_date,
            phase='water_restriction',
            username=username,
            body_weight_g=body_weight_g,
            weight_pct_baseline=weight_pct_baseline,
            weight_delta_from_baseline_g=weight_delta_from_baseline_g,
            bottle_weight_before_g=bottle_weight_before_g,
            bottle_weight_after_g=bottle_weight_after_g,
            water_added_ml=water_added_ml,
            
            task=task,
            task_duration_min=task_duration_min,
            total_trials=total_trials,
            correct_trials=correct_trials,
            correct_rate=correct_rate,
            reward_size_ul=reward_size_ul,
            task_water_ml=task_water_ml,
            remain_water_ml=remain_water_ml,
            
            water_consumed_ml=water_consumed_ml,
            target_water_ml=target_water_ml,
            minimum_water_ml=minimum_water_ml,
            removed_from_restriction=0,
            health_status=health_status,
            notes=notes
        )

    @staticmethod
    def update_previous_consumption(animal_id, current_date, current_bottle_before_g):
        """
        Updates the water_consumed_ml for the previous day's log based on today's bottle retrieval.
        """
        # Find the most recent log before today
        key = {'animal_id': animal_id}
        prev_log = (DailyLog & key & f'log_date < "{current_date}"').fetch(
            'log_date', 'bottle_weight_after_g', 'task_water_ml', 
            order_by='log_date DESC', limit=1, as_dict=True)
            
        if prev_log:
            prev = prev_log[0]
            prev_date = prev['log_date']
            prev_after = float(prev['bottle_weight_after_g']) if prev['bottle_weight_after_g'] else 0.0
            prev_task = float(prev['task_water_ml']) if prev['task_water_ml'] else 0.0
            
            # Consumption = (YesterdayAfter - TodayBefore) + YesterdayTask
            bottle_consumed = max(0, prev_after - float(current_bottle_before_g))
            total_consumed = bottle_consumed + prev_task
            
            # Update the record (DataJoint pattern: fetch, delete, re-insert)
            restriction = key & {'log_date': prev_date}
            record = (DailyLog & restriction).fetch1()
            (DailyLog & restriction).delete_quick()
            record['water_consumed_ml'] = total_consumed
            DailyLog.insert1(record)
            return True, prev_date, total_consumed
            
        return False, None, 0.0


@schema
class WaterCalibration(dj.Manual):
    definition = """
    setup             : varchar(64)   # e.g. "Box 1", "Rig 2"
    calibration_id    : int           # repeat number (1-based index for this setup)
    pump_time_ms      : int           # duration in ms
    ---
    -> experiment.Person
    calibration_time  : datetime      # timestamp of calibration
    continuous_rate_hz: decimal(5,2)  # pump frequency in Hz
    number_of_pulses  : int           # number of pulses delivered
    water_left_ml     : decimal(5,3)  # measured output Left in ml
    water_right_ml    : decimal(5,3)  # measured output Right in ml
    notes=''          : varchar(1024) # comments
    """
