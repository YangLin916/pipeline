import streamlit as st
import pandas as pd
import datetime
import time
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import datajoint as dj

# --- Database Configuration ---
# Credentials should be loaded from pipeline_config.json (gitignored) or environment variables.
# dj.config.save_local() # Optional: saves to local config file for future runs

# Set page config
st.set_page_config(page_title="Active Sensing Daily Log", page_icon="🐭", layout="wide")

# Title
st.title("🐭 Active Sensing Daily Log")

# --- Database Connection ---
try:
    # Check if connection is alive
    try:
        dj.conn().query("SELECT 1")
    except Exception:
        # If dead or shaky, force reconnect
        dj.conn().close()
        dj.conn().connect()
        
    st.success("Connected to Database")
except Exception as e:
    st.error(f"Could not connect to DataJoint: {e}")
    st.stop()

# Import pipeline modules AFTER verification of connection
# This ensures that create_tables=True in schemas works if the table is missing
from pipeline import active_sense, mice

# --- Tabs ---
tab_log, tab_new, tab_calib = st.tabs(["📝 Daily Log", "➕ New Subject", "💧 Water Calibration"])

# ==========================================
# TAB 1: DAILY LOG
# ==========================================
with tab_log:
    # --- Sidebar: Mouse Selection ---
    st.sidebar.header("Select Subject")

    # Fetch active mice with retry logic for connection stability
    try:
        active_mice = active_sense.MouseInfo.fetch('animal_id')
    except Exception:
        # If fetch fails, try reconnecting once
        dj.conn().close()
        dj.conn().connect()
        active_mice = active_sense.MouseInfo.fetch('animal_id')

    if len(active_mice) == 0:
        st.warning("No mice found in MouseInfo table. Please add a new subject in the 'New Subject' tab.")
        selected_mouse = None
    else:
        selected_mouse = st.sidebar.selectbox("Choose Mouse ID", active_mice)

    if selected_mouse:
        # 1. Fetch Mouse Details
        mouse_info = (active_sense.MouseInfo & f'animal_id="{selected_mouse}"').fetch1()

        st.subheader(f"Subject: {selected_mouse}")
        col1, col2, col3 = st.columns(3)
        
        # Calculate target water
        target_water = float(mouse_info['baseline_weight_g']) * 0.05
        min_water = float(mouse_info['baseline_weight_g']) * 0.025
        
        col1.metric("Baseline Weight", f"{mouse_info['baseline_weight_g']} g")
        col2.metric("Target Water (5%)", f"{target_water:.2f} ml")
        col3.metric("Min Water (2.5%)", f"{min_water:.2f} ml")

        st.markdown("---")

        # 2. Data Entry Form
        st.header("New Log Entry")

        with st.container():
            col_date, col_user = st.columns(2)
            log_date = col_date.date_input("Date", datetime.date.today())
            username = col_user.text_input("Experimenter (Username)", value="yang")

            st.markdown("#### Weight Measurements")
            c1, c2, c3 = st.columns(3)
            body_weight = c1.number_input("Body Weight (g)", min_value=0.0, format="%.2f", step=0.1)
            bottle_before = c2.number_input("Bottle Before (Retrieval) (g)", min_value=0.0, format="%.2f", step=0.1)
            bottle_after = c3.number_input("Bottle After (Refill) (g)", min_value=0.0, format="%.2f", step=0.1)

            st.markdown("#### Task Performance")
            
            # Task Selection
            existing_tasks = ["None"]
            try:
                # Fetch distinct tasks (ignore None)
                tasks_db = (active_sense.DailyLog - "task is null").fetch("task")
                if len(tasks_db) > 0:
                     existing_tasks += sorted(list(set(tasks_db)))
            except:
                pass # Table might be empty
            
            # Trigger rerun on change so we can show conditional inputs immediately
            task_select = st.selectbox("Select Task", existing_tasks + ["Enter New Task..."])
            
            task_name = None
            if task_select == "Enter New Task...":
                 task_name = st.text_input("Enter Task Name")
            elif task_select != "None":
                 task_name = task_select
            
            # Task Details (Conditional)
            task_duration = 0
            total_trials = 0
            correct_trials = 0
            reward_size = 0.0
            calc_task_water = 0.0
            
            if task_name:
                t1, t2, t3, t4 = st.columns(4)
                task_duration = t1.number_input("Duration (min)", min_value=0, step=1)
                total_trials = t2.number_input("Total Trials", min_value=0, step=1)
                correct_trials = t3.number_input("Correct Trials", min_value=0, step=1)
                reward_size = t4.number_input("Reward Size (ul)", min_value=0.0, format="%.1f", step=0.5)
                
                # Auto-calc display
                if correct_trials > 0 and reward_size > 0:
                    calc_task_water = (correct_trials * reward_size) / 1000.0
                    st.info(f"💧 Calculated Task Water: **{calc_task_water:.2f} ml** ({correct_trials} * {reward_size}ul)")
                    
                    st.info(f"🚰 Remaining Water Needed: **{max(0, target_water - calc_task_water):.2f} ml**")

            st.markdown("#### Water Management")
            c4, c5 = st.columns(2)
            # Renamed from supplemental -> task_water input (user can override calculated)
            # Since we removed the form, this will update on rerun if calc_task_water changes
            # BUT number_input value=... only sets the initial value.
            # We want it to update IF the user hasn't manually edited it?
            # Streamlit default behavior for value arg is "initial value". 
            # To make it dynamic based on another input, we can use key/state but it's complex.
            # Simpler: Just prompt user to verify. OR we set `value` but it might be ignored if widget exists.
            # Actually, without key change, value is ignored after init. 
            # Let's adding a helper key based on trials to force update? No that resets the widget focus.
            # We will just warn user to check it.
            
            task_water_input = c4.number_input("Task Water (ml)", 
                                             min_value=0.0, value=float(calc_task_water), format="%.2f", step=0.1,
                                             help="Water obtained during behavior task. Defaults to Correct * Reward / 1000.")
            
            water_added = c5.number_input("Water Added (ml) (Manual Top-up)", min_value=0.0, format="%.2f", step=0.1)

            st.markdown("#### Status & Notes")
            health_status = st.select_slider("Health Status (1=Critical, 5=Normal)", options=['1', '2', '3', '4', '5'], value='5')
            notes = st.text_area("Notes")

            submitted = st.button("💾 Save Entry")

            if submitted:
                try:
                    final_task_water = task_water_input
                    # If user left it 0, but we have data, use calculated.
                    if final_task_water == 0 and calc_task_water > 0:
                        final_task_water = calc_task_water
                    
                    # Use the helper method to compute everything
                    entry = active_sense.DailyLog.compute_daily_log(
                        animal_id=selected_mouse,
                        log_date=log_date,
                        body_weight_g=body_weight,
                        bottle_weight_before_g=bottle_before,
                        bottle_weight_after_g=bottle_after,
                        water_added_ml=water_added,
                        
                        task=task_name,
                        task_duration_min=task_duration,
                        total_trials=total_trials,
                        correct_trials=correct_trials,
                        reward_size_ul=reward_size,
                        task_water_ml=final_task_water,
                        
                        username=username,
                        health_status=health_status,
                        notes=notes
                    )
                    
                    # Insert into database
                    active_sense.DailyLog.insert1(entry)
                    st.success(f"Successfully logged entry for {selected_mouse} on {log_date}!")
                    
                    # Update previous day's consumption using today's bottle retrieval
                    updated, prev_date, total_cons = active_sense.DailyLog.update_previous_consumption(
                        selected_mouse, log_date, bottle_before
                    )
                    if updated:
                        st.info(f"🔄 Updated water consumption for {prev_date}: **{total_cons:.2f} ml**")
                        
                    st.balloons()
                    
                except Exception as e:
                    st.error(f"Error saving entry: {e}")

        st.markdown("---")

        # 3. History View
        st.header("📅 Recent History")
        history = (active_sense.DailyLog & f'animal_id="{selected_mouse}"').fetch(
            format="frame", order_by="log_date DESC", limit=7
        ).reset_index()

        if not history.empty:
            # Clean up dataframe for display
            display_cols = [
                'log_date',
                'body_weight_g', 'weight_pct_baseline', 
                'task', 'correct_rate', 'task_water_ml', 'remain_water_ml',
                'water_consumed_ml',
                'health_status'
            ]
            # Define safe formatter
            def safe_fmt(x):
                try:
                    return "{:.2f}".format(float(x))
                except (ValueError, TypeError):
                    return "-"

            st.dataframe(history[display_cols].style.format(safe_fmt, subset=['body_weight_g', 'weight_pct_baseline', 'task_water_ml', 'remain_water_ml', 'water_consumed_ml']))
            
            st.markdown("---")
            with st.expander("🗑️ Manage History (Delete Entries)"):
                st.warning("⚠️ Deleting an entry is permanent.")
                # Dropdown to select date to delete
                # Convert dates to string for selectbox
                date_options = history['log_date'].astype(str).tolist()
                date_to_delete = st.selectbox("Select Date to Delete", options=date_options)
                
                # UI Confirmation to replace terminal prompt
                confirm_delete = st.checkbox(f"✅ I confirm I want to delete the log for {date_to_delete}", key="confirm_delete_check")
                
                if st.button("Delete Selected Entry", type="primary", disabled=not confirm_delete):
                    if date_to_delete and confirm_delete:
                        try:
                            # Toggle safemode in config since delete() arg is not supported
                            # Save current state
                            prev_safemode = dj.config.get('safemode', True)
                            dj.config['safemode'] = False
                            
                            (active_sense.DailyLog & {'animal_id': selected_mouse, 'log_date': date_to_delete}).delete()
                            
                            # Restore safemode
                            dj.config['safemode'] = prev_safemode
                            
                            st.success(f"Deleted entry for {date_to_delete}.")
                            st.rerun() # Refresh to show updated history
                        except Exception as e:
                            st.error(f"Error deleting: {e}")
        else:
            st.info("No logs found for this mouse.")

# ==========================================
# TAB 2: NEW SUBJECT
# ==========================================
# ==========================================
# TAB 2: NEW SUBJECT
# ==========================================
with tab_new:
    st.header("Add New Subject")
    st.info("Enter an Animal ID to check if it exists in the central `mice` database.")

    # Input ID outside the form to allow dynamic lookup
    check_animal_id = st.text_input("Animal ID", key="new_subject_id")
    
    found_mouse = None
    if check_animal_id:
        # Check central database
        found_mouse = (mice.Mice & f'animal_id="{check_animal_id}"').fetch1() if (mice.Mice & f'animal_id="{check_animal_id}"') else None
    
    if found_mouse:
        st.success(f"✅ Found **{check_animal_id}** in central database. Pre-filling known info.")
        known_sex = found_mouse.get('sex', 'U')
        known_dob = found_mouse.get('dob', datetime.date.today())
        
        # When found, we skip generating a "create mouse" step and go straight to enrollment with manual extras
        with st.form("enroll_existing_form", enter_to_submit=False):
            st.write(f"**Sex:** {known_sex}")
            st.write(f"**DOB:** {known_dob}")
            
            st.subheader("Additional Active Sense Info")
            c1, c2 = st.columns(2)
            line = c1.text_input("Line / Strain (Required)", value="C57BL/6J")
            genotype = c2.text_input("Genotype", value="")
            
            name_input = st.text_input("Name (Optional)", value="", help="Pet name or alias")
            
            c3, c4 = st.columns(2)
            single_housing = c3.date_input("Single Housing Date", datetime.date.today())
            restriction_start = c4.date_input("Restriction Start Date", datetime.date.today() + datetime.timedelta(days=3))
            
            base_weight = st.number_input("Baseline Weight (g)", min_value=0.0, format="%.2f", step=0.1)
            notes = st.text_area("Notes")

            submit_enroll = st.form_submit_button("Enroll Subject")

            if submit_enroll:
                try:
                    active_sense.MouseInfo.insert1(dict(
                        animal_id=check_animal_id,
                        sex=known_sex,
                        dob=known_dob,
                        name=name_input,
                        line=line,
                        genotype=genotype,
                        single_housing_date=single_housing,
                        restriction_start_date=restriction_start,
                        baseline_weight_g=base_weight,
                        baseline_weight_date=datetime.date.today(),
                        notes=notes
                    ), skip_duplicates=True)
                    st.success(f"Successfully enrolled {check_animal_id} ({name_input})!")
                    st.info("Refresh to see in Daily Log.")
                except Exception as e:
                    st.error(f"Error enrolling: {e}")

    elif check_animal_id:
        st.warning(f"⚠️ **{check_animal_id}** not found. Please enter FULL details to create it.")
        
        with st.form("create_full_form", enter_to_submit=False):
            st.subheader("Basic Info")
            c1, c2, c3 = st.columns(3)
            sex = c1.selectbox("Sex", ["M", "F", "U"])
            dob = c2.date_input("Date of Birth", datetime.date.today())
            color = c3.selectbox("Color", ["Black", "Brown", "White", "unknown"])
            
            c4, c5 = st.columns(2)
            ear_punch = c4.selectbox("Ear Punch", ["None", "R", "L", "RL", "RR", "LL", "unknown"])
            line = c5.text_input("Line / Strain", value="C57BL/6J")
            
            st.subheader("Active Sense Info")
            c6, c7 = st.columns(2)
            genotype = c6.text_input("Genotype", value="")
            name_input = c7.text_input("Name (Optional)", value="")

            base_weight = st.number_input("Baseline Weight (g)", min_value=0.0, format="%.2f", step=0.1)
            
            c8, c9 = st.columns(2)
            single_housing = c8.date_input("Single Housing Date", datetime.date.today())
            restriction_start = c9.date_input("Restriction Start Date", datetime.date.today() + datetime.timedelta(days=3))
            
            notes = st.text_area("Initial Notes")

            submit_new = st.form_submit_button("➕ Create & Enroll")

            if submit_new:
                try:
                    # 1. Insert into central mice
                    mice.Mice.insert1(dict(
                        animal_id=check_animal_id,
                        sex=sex,
                        dob=dob,
                        color=color,
                        ear_punch=ear_punch
                    ), skip_duplicates=True)
                    
                    # 2. Insert into MouseInfo
                    # Try to insert 'name' separately/gracefully or assume schema updated
                    active_sense.MouseInfo.insert1(dict(
                        animal_id=check_animal_id,
                        sex=sex,
                        dob=dob,
                        name=name_input,
                        line=line,
                        genotype=genotype,
                        single_housing_date=single_housing,
                        restriction_start_date=restriction_start,
                        baseline_weight_g=base_weight,
                        baseline_weight_date=datetime.date.today(),
                        notes=notes
                    ), skip_duplicates=True)
                    
                    st.success(f"Successfully created and enrolled {check_animal_id}!")
                    st.info("Refresh to see in Daily Log.")
                except Exception as e:
                    st.error(f"Error creating: {e}")

# ==========================================
# TAB 3: WATER CALIBRATION
# ==========================================
with tab_calib:
    st.header("Water Calibration Log")
    
    # 1. Input Form
    with st.expander("➕ Add New Calibration", expanded=True):
        col1, col2 = st.columns(2)
        setup_name = col1.text_input("Setup Name", value="Box 1", help="e.g. Box 1, Rig 2")
        calib_user = col2.text_input("Experimenter", value="yang")
        
        # Auto-calculate Repeat ID based on Setup
        repeat_id = 1
        if setup_name:
            try:
                # Find max calibration_id for this setup
                existing_ids = (active_sense.WaterCalibration & f'setup="{setup_name}"').fetch('calibration_id')
                if len(existing_ids) > 0:
                    repeat_id = max(existing_ids) + 1
            except:
                pass
        
        st.info(f"🔢 Next Calibration ID for **{setup_name}**: **{repeat_id}**")
        
        with st.form("water_calib_form", enter_to_submit=False):
            st.subheader("Pump Parameters")
            c1, c2, c3 = st.columns(3)
            pump_time = c1.number_input("Pump Time (ms)", min_value=0, step=10, value=150)
            cont_rate = c2.number_input("Continuous Rate (Hz)", min_value=0.0, step=1.0, value=20.0, format="%.2f")
            pulses = c3.number_input("Number of Pulses", min_value=1, step=1, value=100)
            
            st.subheader("Measurement")
            note_col, water_c1, water_c2 = st.columns([2, 1, 1])
            notes = note_col.text_area("Notes", placeholder="e.g. 500 pulses total, measured with cylinder...")
            water_left = water_c1.number_input("Left Volume (ml)", min_value=0.0, format="%.3f", step=0.05)
            water_right = water_c2.number_input("Right Volume (ml)", min_value=0.0, format="%.3f", step=0.05)
            
            submit_calib = st.form_submit_button("💾 Save Calibration")
            
            if submit_calib:
                try:
                    active_sense.WaterCalibration.insert1(dict(
                        setup=setup_name,
                        calibration_id=repeat_id,
                        username=calib_user,
                        calibration_time=datetime.datetime.now(),
                        pump_time_ms=pump_time,
                        continuous_rate_hz=cont_rate,
                        number_of_pulses=pulses,
                        water_left_ml=water_left,
                        water_right_ml=water_right,
                        notes=notes
                    ))
                    st.success(f"Saved calibration #{repeat_id} for {setup_name}!")
                    st.rerun()
                except Exception as e:
                    st.error(f"Error saving calibration: {e}")

    # 2. History
    st.markdown("---")
    st.subheader("📜 Calibration History")
    
    # Filter by setup
    all_setups = ["All"]
    try:
        db_setups = (active_sense.WaterCalibration).fetch('setup')
        if len(db_setups) > 0:
            all_setups += sorted(list(set(db_setups)))
    except:
        pass
        
    filter_setup = st.selectbox("Filter by Setup", all_setups)
    
    try:
        if filter_setup == "All":
            calib_history = active_sense.WaterCalibration.fetch(format="frame", order_by="calibration_time DESC", limit=20).reset_index()
        else:
            calib_history = (active_sense.WaterCalibration & f'setup="{filter_setup}"').fetch(format="frame", order_by="calibration_time DESC", limit=20).reset_index()
            
        if not calib_history.empty:
            st.dataframe(calib_history[['calibration_time', 'setup', 'calibration_id', 'username', 'pump_time_ms', 'water_left_ml', 'water_right_ml', 'notes']])
        else:
            st.info("No calibration records found.")
    except Exception as e:
        st.warning(f"Could not fetch history (Table might be empty or missing): {e}")

    # 3. Analysis & Calculator
    st.markdown("---")
    with st.expander("📈 Calibration Analysis & Calculator", expanded=False):
        if filter_setup != "All":
            # Fetch data for this setup
            try:
                df_calib = (active_sense.WaterCalibration & f'setup="{filter_setup}"').fetch(format="frame").reset_index()
                
                if len(df_calib) > 2:
                    # Function to fit: Time = m * Volume + c
                    def func_linear(x, m, c):
                        return m * x + c

                    x_left = df_calib['water_left_ml'].values.astype(float)
                    x_right = df_calib['water_right_ml'].values.astype(float)
                    y_time = df_calib['pump_time_ms'].values.astype(float)

                    # Fit Left
                    try:
                        popt_l, _ = curve_fit(func_linear, x_left, y_time)
                        m_l, c_l = popt_l
                        valid_l = True
                    except:
                        valid_l = False
                    
                    # Fit Right
                    try:
                        popt_r, _ = curve_fit(func_linear, x_right, y_time)
                        m_r, c_r = popt_r
                        valid_r = True
                    except:
                        valid_r = False

                    # Calculator UI
                    st.subheader("💧 Reward Calculator")
                    target_vol_ul = st.number_input("Target Reward Size (µl)", value=5.0, step=0.5)
                    target_vol_ml = target_vol_ul / 1000.0
                    
                    col_res1, col_res2 = st.columns(2)
                    
                    if valid_l:
                        req_time_l = func_linear(target_vol_ml, *popt_l)
                        col_res1.success(f"**Left Port**: {req_time_l:.1f} ms")
                        col_res1.caption(f"Fit: T = {m_l:.1f}*V + {c_l:.1f}")
                    else:
                        col_res1.warning("Left: Not enough data to fit")

                    if valid_r:
                        req_time_r = func_linear(target_vol_ml, *popt_r)
                        col_res2.success(f"**Right Port**: {req_time_r:.1f} ms")
                        col_res2.caption(f"Fit: T = {m_r:.1f}*V + {c_r:.1f}")
                    else:
                        col_res2.warning("Right: Not enough data to fit")

                    # plotting
                    st.subheader("Curve Fit")
                    fig, ax = plt.subplots(figsize=(6, 4))
                    
                    # Plot raw data
                    ax.scatter(x_left, y_time, color='blue', label='Left Data', alpha=0.6)
                    ax.scatter(x_right, y_time, color='red', label='Right Data', alpha=0.6)
                    
                    # Plot fits
                    x_plot = np.linspace(0, max(x_left.max(), x_right.max()) * 1.1, 50)
                    if valid_l:
                        ax.plot(x_plot, func_linear(x_plot, *popt_l), 'b--', alpha=0.5, label='Left Fit')
                    if valid_r:
                        ax.plot(x_plot, func_linear(x_plot, *popt_r), 'r--', alpha=0.5, label='Right Fit')
                        
                    ax.set_xlabel('Water Volume (ml)')
                    ax.set_ylabel('Pump Time (ms)')
                    ax.legend()
                    ax.grid(True, linestyle=':', alpha=0.6)
                    
                    st.pyplot(fig)
                    plt.close(fig) # Close to avoid memory leaks
                    
                else:
                    st.warning("Need at least 3 data points for this setup to perform analysis.")
            except Exception as e:
                st.error(f"Analysis Error: {e}")
        else:
            st.info("Please select a specific 'Setup' filter above to enable analysis.")
