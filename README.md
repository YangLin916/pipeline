# ninai/pipeline
Processing pipeline for scans and behavioral files.

# Docker container
* [ninai/pipeline](https://hub.docker.com/r/ninai/pipeline/) `ninai/pipeline`

# Local Installation & Active Sense App

## Prerequisites
- Python 3.8+
- [pip](https://pip.pypa.io/en/stable/installation/)
- [git](https://git-scm.com/downloads)

## Installation on a New Computer

1. **Clone the repository:**
   ```bash
   git clone <repository_url>
   cd pipeline
   ```

2. **Create and activate a virtual environment (Recommended):**
   ```bash
   python3 -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies:**
   You can install from `requirements.txt` or install the package in editable mode.
   ```bash
   pip install -r requirements.txt
   # OR
   pip install -e python/
   ```
   *Note: This pipeline depends on the `commons` package from the atlab organization. `requirements.txt` handles this automatically via git.*

4. **Configuration:**
   Copy the example configuration file and fill in your database credentials.
   ```bash
   cp python/pipeline_config.example.json python/pipeline_config.json
   ```
   **Edit `python/pipeline_config.json`:**
   - Update `database.host`, `database.user`, and `database.password` with your DataJoint credentials.
   - Adjust `path.mounts` if necessary.

## Running the Active Sense Daily Log App

To launch the daily log application:

```bash
streamlit run python/apps/daily_log_app.py
```

The app should open automatically in your browser.
