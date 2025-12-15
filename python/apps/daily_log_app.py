import streamlit as st
import pandas as pd
import datetime
import datajoint as dj

# --- Database Configuration ---
# Credentials should be loaded from pipeline_config.json (gitignored) or environment variables.
# dj.config.save_local() # Optional: saves to local config file for future runs

from pipeline import active_sense, mice

# Set page config
st.set_page_config(page_title="Active Sense Daily Log", page_icon="🐭", layout="wide")

# Title
st.title("🐭 Active Sense Daily Restriction Log")

# --- Database Connection ---
try:
    dj.conn()
except Exception as e:
    st.error(f"Could not connect to DataJoint: {e}")
    st.stop()

# --- Tabs ---
tab_log, tab_new = st.tabs(["📝 Daily Log", "➕ New Subject"])

# ==========================================
# TAB 1: DAILY LOG
# ==========================================
with tab_log:
    # --- Sidebar: Mouse Selection ---
    st.sidebar.header("Select Subject")

    # Fetch active mice
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
        col1.metric("Baseline Weight", f"{mouse_info['baseline_weight_g']} g")
        col2.metric("Target Water (5%)", f"{float(mouse_info['baseline_weight_g']) * 0.05:.1f} ml")
        col3.metric("Min Water (2.5%)", f"{float(mouse_info['baseline_weight_g']) * 0.025:.1f} ml")

        st.markdown("---")

        # 2. Data Entry Form
        st.header("New Log Entry")

        with st.form("daily_log_form"):
            col_date, col_user = st.columns(2)
            log_date = col_date.date_input("Date", datetime.date.today())
            username = col_user.text_input("Experimenter (Username)", value="yang")

            st.markdown("#### Weight Measurements")
            c1, c2, c3 = st.columns(3)
            body_weight = c1.number_input("Body Weight (g)", min_value=0.0, format="%.2f", step=0.1)
            bottle_before = c2.number_input("Bottle Before (Retrieval) (g)", min_value=0.0, format="%.2f", step=0.1)
            bottle_after = c3.number_input("Bottle After (Refill) (g)", min_value=0.0, format="%.2f", step=0.1)

            st.markdown("#### Water Management")
            c4, c5 = st.columns(2)
            water_added = c4.number_input("Water Added (ml)", min_value=0.0, format="%.2f", step=0.1)
            supp_water = c5.number_input("Supplemental Water (ml)", min_value=0.0, value=0.0, format="%.2f", step=0.1)

            st.markdown("#### Status & Notes")
            health_status = st.select_slider("Health Status (1=Critical, 5=Normal)", options=['1', '2', '3', '4', '5'], value='5')
            notes = st.text_area("Notes")

            submitted = st.form_submit_button("💾 Save Entry")

            if submitted:
                try:
                    # Use the helper method to compute everything
                    entry = active_sense.DailyRestrictionLog.compute_daily_log(
                        animal_id=selected_mouse,
                        log_date=log_date,
                        body_weight_g=body_weight,
                        bottle_weight_before_g=bottle_before,
                        bottle_weight_after_g=bottle_after,
                        water_added_ml=water_added,
                        supplemental_water_ml=supp_water,
                        username=username,
                        health_status=health_status,
                        notes=notes
                    )
                    
                    # Insert into database
                    active_sense.DailyRestrictionLog.insert1(entry)
                    st.success(f"Successfully logged entry for {selected_mouse} on {log_date}!")
                    st.balloons()
                    
                except Exception as e:
                    st.error(f"Error saving entry: {e}")

        st.markdown("---")

        # 3. History View
        st.header("📅 Recent History")
        history = (active_sense.DailyRestrictionLog & f'animal_id="{selected_mouse}"').fetch(
            format="frame", order_by="log_date DESC", limit=7
        )

        if not history.empty:
            # Clean up dataframe for display
            display_cols = [
                'body_weight_g', 'weight_pct_baseline', 
                'water_consumed_ml', 'bottle_weight_before_g', 'bottle_weight_after_g',
                'health_status', 'notes'
            ]
            st.dataframe(history[display_cols].style.format("{:.2f}", subset=['body_weight_g', 'weight_pct_baseline', 'water_consumed_ml', 'bottle_weight_before_g', 'bottle_weight_after_g']))
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
        with st.form("enroll_existing_form"):
            st.write(f"**Sex:** {known_sex}")
            st.write(f"**DOB:** {known_dob}")
            
            st.subheader("Additional Active Sense Info")
            c1, c2 = st.columns(2)
            line = c1.text_input("Line / Strain (Required)", value="C57BL/6J")
            genotype = c2.text_input("Genotype", value="")
            
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
                        line=line,
                        genotype=genotype,
                        single_housing_date=single_housing,
                        restriction_start_date=restriction_start,
                        baseline_weight_g=base_weight,
                        baseline_weight_date=datetime.date.today(),
                        notes=notes
                    ), skip_duplicates=True)
                    st.success(f"Successfully enrolled {check_animal_id}!")
                    st.info("Refresh to see in Daily Log.")
                except Exception as e:
                    st.error(f"Error enrolling: {e}")

    elif check_animal_id:
        st.warning(f"⚠️ **{check_animal_id}** not found. Please enter FULL details to create it.")
        
        with st.form("create_full_form"):
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
            base_weight = c7.number_input("Baseline Weight (g)", min_value=0.0, format="%.2f", step=0.1)
            
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
                    active_sense.MouseInfo.insert1(dict(
                        animal_id=check_animal_id,
                        sex=sex,
                        dob=dob,
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
