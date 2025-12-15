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
    line=''                               : varchar(128)              # strain or line name
    genotype=''                           : varchar(255)              # genotype details
    single_housing_date = null            : date                      # start of single housing (baseline days 1-3)
    restriction_start_date = null         : date                      # start of water restriction (days 4+)
    baseline_weight_g = null              : decimal(6,3)              # average baseline weight (g)
    baseline_weight_date = null           : date                      # date baseline weight was measured
    notes=''                              : varchar(1024)             # free-text notes
    """


@schema
class DailyRestrictionLog(dj.Manual):
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
    water_added_ml = null               : decimal(6,3)                  # amount of water added (ml), including dish deliveries
    supplemental_water_ml = 0           : decimal(6,3)                  # extra water given (e.g. syringe/gavage)
    water_consumed_ml = null            : decimal(6,3)                  # calculated intake (bottle delta + supplemental)
    target_water_ml = null              : decimal(6,3)                  # allowance based on baseline weight * 50 ml/kg
    minimum_water_ml = null             : decimal(6,3)                  # minimum allowance based on 25 ml/kg rule
    removed_from_restriction = 0        : boolean                       # mark if mouse taken off restriction today
    health_status = '5'                 : enum('1', '2', '3', '4', '5') # health score (1=Critical, 5=Normal)
    notes = ''                          : varchar(1024)                 # free-text notes (spillage, handling, etc.)
    """

    @staticmethod
    def compute_daily_log(animal_id, log_date, body_weight_g,
                          bottle_weight_before_g, bottle_weight_after_g,
                          water_added_ml, supplemental_water_ml=0,
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

        # Calculate consumption (requires previous day's bottle_weight_after)
        # We look for the most recent log before this date
        prev_logs = (DailyRestrictionLog & {'animal_id': animal_id} & f'log_date < "{log_date}"').fetch(
            'bottle_weight_after_g', order_by='log_date DESC', limit=1)
        
        water_consumed_ml = None
        if len(prev_logs) > 0:
            prev_after = float(prev_logs[0])
            water_consumed_ml = (prev_after - bottle_weight_before_g) + supplemental_water_ml

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
            supplemental_water_ml=supplemental_water_ml,
            water_consumed_ml=water_consumed_ml,  # calculated
            target_water_ml=target_water_ml,
            minimum_water_ml=minimum_water_ml,
            removed_from_restriction=0,
            health_status=health_status,
            notes=notes
        )
