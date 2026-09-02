"""
CalTop TAU Profiling Time-Share Visualization
=============================================

This script parses a TAU profile generated from a CalTop execution and
constructs a pie chart showing the relative computational cost of the
major stages of the topology optimization workflow.

TAU inclusive timer values are extracted directly from the profile file
and converted from microseconds to seconds. The measured execution time
is grouped into the following categories:

    - Linear solve:
        PARDISO time associated with the primal structural analysis.

    - Sensitivity:
        Inclusive time associated with the stress, compliance,
        center-of-gravity (CG), and volume sensitivity evaluations.
        The stress sensitivity timer includes the PARDISO solve required
        for the stress-adjoint system.

    - Preprocessing:
        CalTop preprocessing and linear-system assembly.

    - Filter:
        Time associated with the CalTop filtering operations.

    - Output:
        Time associated with CalTop file I/O.

The current TAU instrumentation reports two calls to the
"PARDISO: Total" timer: one for the primal structural system and one for
the stress-adjoint system. Because TAU reports their aggregate inclusive
time rather than the individual duration of each call, this script
temporarily assumes that the two calls contribute equally to the total
PARDISO time. Consequently, one half of the aggregate PARDISO time is
assigned to the Linear solve category. The second PARDISO call is not
added separately to the Sensitivity category because it is already
contained within the inclusive "Stress Grad: Total" timer.

The script prints the timing breakdown and corresponding percentages to
the terminal and generates a publication-quality pie chart using the
Seaborn colorblind palette and Times New Roman font.

Input
-----
TAU profile file:
    PROFILE_FILE, typically of the form "profile.0.0.0".

Output
------
Figure:
    CalTop_time_share.png

Notes
-----
All timing quantities used for the computational breakdown are TAU
inclusive times. Because inclusive timers may contain nested timers,
care must be taken to avoid double counting when constructing the
aggregate timing categories.
"""


import re
import matplotlib.pyplot as plt
import seaborn as sns


# ============================================================
# User settings
# ============================================================

PROFILE_FILE = "Threads_1/profile.0.0.0"

FONT_NAME = "Times New Roman"
FONT_SIZE = 22
PERCENT_FONT_SIZE = 26

FIGSIZE = (8, 6)
DPI = 300

OUTPUT_FILE = "CalTop_time_share.png"


# ============================================================
# Matplotlib / Seaborn settings
# ============================================================

sns.set_palette("colorblind")

plt.rcParams["font.family"] = FONT_NAME
plt.rcParams["font.size"] = 24
plt.rcParams["axes.labelsize"] = FONT_SIZE
plt.rcParams["xtick.labelsize"] = FONT_SIZE
plt.rcParams["ytick.labelsize"] = FONT_SIZE
plt.rcParams["legend.fontsize"] = 26


# ============================================================
# Read TAU profile
# ============================================================

def read_tau_profile(profile_file):
    """
    Read a TAU profile file.

    Returns:
        times : dictionary containing inclusive time [s]
        calls : dictionary containing number of calls

    TAU profile columns:
        Name  Calls  Subrs  Excl  Incl  ProfileCalls

    Excl and Incl are stored in microseconds.
    """

    times = {}
    calls = {}

    with open(profile_file, "r") as f:

        for line in f:

            line = line.strip()

            if not line.startswith('"'):
                continue

            match = re.match(
                r'^"([^"]+)"\s+'
                r'(\d+)\s+'
                r'(\d+)\s+'
                r'([\d.eE+-]+)\s+'
                r'([\d.eE+-]+)\s+'
                r'(\d+)',
                line
            )

            if match is None:
                continue

            name = match.group(1)

            number_calls = int(match.group(2))

            inclusive_us = float(match.group(5))

            # Convert microseconds -> seconds
            inclusive_s = inclusive_us * 1.0e-6

            times[name] = inclusive_s
            calls[name] = number_calls

    return times, calls


# ============================================================
# Read profile
# ============================================================

times, calls = read_tau_profile(PROFILE_FILE)


# ============================================================
# Helper functions
# ============================================================

def get_timer(name):
    """
    Return inclusive timer value.
    """

    if name not in times:

        raise KeyError(
            f'Timer "{name}" was not found in {PROFILE_FILE}'
        )

    return times[name]


def get_calls(name):
    """
    Return number of calls for a timer.
    """

    if name not in calls:

        raise KeyError(
            f'Timer "{name}" was not found in {PROFILE_FILE}'
        )

    return calls[name]


# ============================================================
# Individual TAU inclusive times
# ============================================================

preprocess = get_timer(
    "preprocess: CalTop"
)

assembly = get_timer(
    "Linear system: Assembly"
)

pardiso_total = get_timer(
    "PARDISO: Total"
)

pardiso_calls = get_calls(
    "PARDISO: Total"
)

file_io = get_timer(
    "fileIO: CalTop"
)

filter_caltop = get_timer(
    "filter: CalTop"
)

stress_grad = get_timer(
    "Stress Grad: Total"
)

compliance_grad = get_timer(
    "Compliance Grad: Total"
)

cg_grad = get_timer(
    "CG Grad: Total"
)

vol_grad = get_timer(
    "Vol Grad: Total"
)


# ============================================================
# Split PARDISO time
# ============================================================

#
# Current TAU profile contains two PARDISO calls:
#
#     Call 1 : primal structural solve
#     Call 2 : stress-adjoint solve
#
# TAU reports only the aggregate PARDISO: Total inclusive time.
#
# Therefore, until the two calls are instrumented separately,
# divide the total time equally between the two calls.
#

if pardiso_calls != 2:

    raise ValueError(
        "\nExpected PARDISO: Total to have exactly 2 calls, "
        f"but found {pardiso_calls}.\n"
        "The current PARDISO splitting assumption is therefore "
        "not valid."
    )


pardiso_primal = pardiso_total / 2.0

pardiso_stress_adjoint = pardiso_total / 2.0


# ============================================================
# Construct plotting categories
# ============================================================

# ------------------------------------------------------------
# Preprocessing
#
# CalTop preprocessing + linear-system assembly
# ------------------------------------------------------------

preprocessing_total = (
    preprocess
    + assembly
)


# ------------------------------------------------------------
# Linear solve
#
# Only the FIRST PARDISO call is assigned here.
#
# The second PARDISO call is the stress-adjoint solve.
# ------------------------------------------------------------

linear_solve_total = pardiso_primal


# ------------------------------------------------------------
# Sensitivity
#
# Stress Grad: Total is an INCLUSIVE timer.
#
# It already contains the second PARDISO call used for the
# stress adjoint. Therefore, we do NOT add
# pardiso_stress_adjoint again here.
# ------------------------------------------------------------

sensitivity_total = (
    stress_grad
    + compliance_grad
    + cg_grad
    + vol_grad
)


# ------------------------------------------------------------
# Filter
#
# Total CalTop filtering stage
# ------------------------------------------------------------

filter_total = filter_caltop


# ------------------------------------------------------------
# Output
#
# CalTop file I/O
# ------------------------------------------------------------

output_total = file_io


# ============================================================
# Plotting data
# ============================================================

labels = [
    "Linear solve",
    "Sensitivity",
    "Preprocessing",
    "Filter",
    "Output",
]

values = [
    linear_solve_total,
    sensitivity_total,
    preprocessing_total,
    filter_total,
    output_total,
]


# ============================================================
# Print detailed timing information
# ============================================================

print()
print("PARDISO timing")
print("========================================")

print(
    f"PARDISO total       : "
    f"{pardiso_total:10.6f} s"
)

print(
    f"PARDISO calls       : "
    f"{pardiso_calls:d}"
)

print(
    f"Primal solve        : "
    f"{pardiso_primal:10.6f} s"
)

print(
    f"Stress adjoint      : "
    f"{pardiso_stress_adjoint:10.6f} s"
)

print()
print(
    "NOTE: primal and adjoint PARDISO times are currently "
    "estimated by dividing the aggregate PARDISO time equally."
)


# ============================================================
# Print pie-chart timing summary
# ============================================================

total = sum(values)

print()
print("Inclusive timing summary")
print("========================================")

for label, value in zip(labels, values):

    percentage = 100.0 * value / total

    print(
        f"{label:<15s}: "
        f"{value:10.6f} s "
        f"({percentage:6.2f} %)"
    )

print("----------------------------------------")

print(
    f"{'Total':<15s}: "
    f"{total:10.6f} s"
)


# ============================================================
# Pie chart
# ============================================================

fig, ax = plt.subplots(
    figsize=FIGSIZE
)


# ============================================================
# Seaborn colorblind palette
# ============================================================

colors = sns.color_palette(
    "colorblind",
    n_colors=len(values)
)


# ============================================================
# Create pie chart
# ============================================================

wedges, texts, autotexts = ax.pie(

    values,

    # No labels directly on pie
    labels=None,

    # Seaborn colorblind palette
    colors=colors,

    # Percentage format
    autopct="%1.2f%%",

    # Percentage labels inside pie
    pctdistance=0.70,

    # Pie orientation
    startangle=90,
    counterclock=False,

    # Pie radius
    radius=0.90,

    # Wedge appearance
    wedgeprops={
        "edgecolor": "white",
        "linewidth": 1.5,
        "alpha": 1.0
    },

    # Percentage text properties
    textprops={
        "fontsize": PERCENT_FONT_SIZE,
        "fontfamily": FONT_NAME
    }
)


# ============================================================
# Percentage labels
# ============================================================

for text in autotexts:

    text.set_fontfamily(FONT_NAME)
    text.set_fontsize(PERCENT_FONT_SIZE)

    text.set_color("black")


# ============================================================
# Legend
#
# Right side of pie, single column
# ============================================================

ax.legend(

    wedges,
    labels,

    # Place legend to right of pie
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),

    # Single column
    ncols=1,

    frameon=True,

    # Legend spacing
    labelspacing=1.0,
    handletextpad=0.5,

    prop={
        "family": FONT_NAME,
        "size": 24
    }
)


# ============================================================
# Keep pie circular
# ============================================================

ax.axis("equal")


# ============================================================
# Figure size
# ============================================================

F = plt.gcf()

Size = F.get_size_inches()

F.set_size_inches(
    Size[0] * 1.5,
    Size[1] * 1.5,
    forward=True
)


# ============================================================
# Layout
# ============================================================

plt.tight_layout()


# ============================================================
# Save figure
# ============================================================

plt.savefig(
    OUTPUT_FILE,
    dpi=DPI,
    bbox_inches="tight"
)


# ============================================================
# Display
# ============================================================

plt.show()