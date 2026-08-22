"""BioProtein Studio analysis modules.

On Streamlit Community Cloud, a few bioinformatics executables are not
available through Debian apt. We bootstrap only those missing tools into a
user-local cache. Existing local/Conda installations always take priority.
"""

try:
    from .gdm_cloud_tools import bootstrap_cloud_tools
    bootstrap_cloud_tools()
except Exception:
    # Never block the app if cloud bootstrap fails; individual modules expose
    # tool readiness and manual/reference fallbacks.
    pass
